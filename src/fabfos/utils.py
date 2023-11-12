import os
import time
from datetime import datetime as dt
from pathlib import Path
import re

USER = "hallamlab" # github id
MODULE_ROOT = Path("/".join(os.path.realpath(__file__).split('/')[:-1]))
NAME = MODULE_ROOT.name.lower()
ENTRY_POINTS = [NAME, "ffs"]

def _get_version() -> str:
    with open(MODULE_ROOT.joinpath("version.txt")) as v:
        return v.readline()
VERSION = _get_version()

def regex(r, s):
    for m in re.finditer(r, s):
        yield s[m.start():m.end()]

class StdTime:
    FORMAT = '%Y-%m-%d_%H-%M-%S'

    @classmethod
    def Timestamp(cls, timestamp: dt|None = None):
        ts = dt.now() if timestamp is None else timestamp
        return f"{ts.strftime(StdTime.FORMAT)}"
    
    @classmethod
    def Parse(cls, timestamp: str|int):
        if isinstance(timestamp, str):
            return dt.strptime(timestamp, StdTime.FORMAT)
        else:
            return dt.fromtimestamp(timestamp/1000)
    
    @classmethod
    def CurrentTimeMillis(cls):
        return round(time.time() * 1000)
    
