import os
import re
import signal
from dataclasses import dataclass
from typing import IO, Any, Callable
from threading import Condition
from threading import Thread
from multiprocessing import Queue
import subprocess
import select

@dataclass
class ShellResult:
    killed: bool
    exit_code: int|None

# example: colors, escape, control sequences
# https://stackoverflow.com/questions/14693701/how-can-i-remove-the-ansi-escape-sequences-from-a-string-in-python 
def StripANSI(s: str):
    return re.compile(r'\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])').sub('', s)

def Shell(cmd: str, on_out: Callable[[str], Any]|None=None, on_err: Callable[[str], Any]|None=None):
    killed = False
    with LiveShell() as shell:
        if on_out is not None: shell.RegisterOnOut(lambda x: on_out(shell.Decode(x)))
        if on_err is not None: shell.RegisterOnErr(lambda x: on_err(shell.Decode(x)))
        try:
            shell.Write(cmd)
            shell.Write(f"exit")
            shell._console.wait()
        except KeyboardInterrupt:
            killed = True
            try:
                os.killpg(os.getpgid(shell.pid), signal.SIGTERM)
            except ProcessLookupError:
                pass
    return ShellResult(killed, shell._console.poll())

class LiveShell:
    class Pipe:
        def __init__(self, io:IO[bytes]|None, lock: Condition=Condition(), q: Queue=Queue()) -> None:
            assert io is not None
            self.IO = io
            self.Lock = lock
            self.Q = q

        def __enter__(self):
            self.Lock.acquire()

        def __exit__(self, exc_type, exc_val, exc_tb):
            self.Lock.release()

    def __init__(self) -> None:
        import pty

        # https://stackoverflow.com/questions/41542960/run-interactive-bash-with-popen-and-a-dedicated-tty-python
        out_master, out_slave = pty.openpty()
        err_master, err_slave = pty.openpty()
        # self._fds = [out_master, err_master, out_slave, err_slave]
        self._fds = [out_master, err_master]

        console = subprocess.Popen(
            ["/bin/bash"],
            stdin=subprocess.PIPE,
            stdout=out_slave,
            stderr=err_slave,
            close_fds=True
        )

        self.ENCODING = "utf-8"
        self._console = console
        self._in = LiveShell.Pipe(console.stdin)
        self._onCloseLock = Condition()
        self._closed = False
        self.pid = console.pid
        self._on_out_callbacks = []
        self._on_err_callbacks = []

        workers: list[Thread] = []
        def reader(fd: int, callbacks):
            _buffer = []
            def _try_read():
                nonlocal _buffer
                changed = False
                while True:
                    # https://stackoverflow.com/a/21429655/13690762
                    r, _, _ = select.select([ fd ], [], [], 0.1)
                    if fd in r:
                        _buffer.append(os.read(fd, 1024))
                        changed = True
                    else:
                        break
                if not changed: return

                for line in b''.join(_buffer).splitlines(True):
                    if line.endswith(b'\n'):
                        yield line
                    else:
                        _buffer = [line]
                        break

            while True:
                if self.IsClosed(): break
                try:
                    for line in _try_read():
                        if line is None: break
                        for cb in callbacks: cb(line)
                except OSError: # fd closed
                    break

        workers.append(Thread(target=reader, args=[out_master, self._on_out_callbacks]))
        workers.append(Thread(target=reader, args=[err_master, self._on_err_callbacks]))
        self._workers = workers
        for w in workers:
            # w.daemon = True # stop with program
            w.start()

    def Send(self, payload: bytes):
        stdin = self._in
        with self._in:
            stdin.IO.write(payload)
            stdin.IO.flush()
    
    def Decode(self, payload: bytes):
        return payload.decode(encoding=self.ENCODING)

    def Write(self, msg: str):
        self.Send(bytes('%s\n' % (msg), encoding=self.ENCODING))

    def RegisterOnOut(self, callback: Callable[[bytes], None]):
        self._on_out_callbacks.append(callback)

    def RegisterOnErr(self, callback: Callable[[bytes], None]):
        self._on_err_callbacks.append(callback)

    def RemoveOnOut(self, callback: Callable[[bytes], None]):
        self._on_out_callbacks.remove(callback)

    def RemoveOnErr(self, callback: Callable[[bytes], None]):
        self._on_err_callbacks.remove(callback)
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.Dispose()
        return

    def IsClosed(self):
        with self._onCloseLock:
            return self._closed

    def Dispose(self):
        with self._onCloseLock:
            if self._closed:
                return
            self._closed = True
            self._onCloseLock.notify_all()

        for w in self._workers:
            w.join()

        self._console.terminate()
        for i, fd in enumerate(self._fds):
            os.close(fd)