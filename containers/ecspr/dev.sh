NAME=ecspr
DOCKER_IMAGE=quay.io/hallamlab/external_$NAME
VERSION=2026.07.14
echo image: $DOCKER_IMAGE:$VERSION
echo ""

HERE=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

case $1 in
    --build|-b)
        # pre-download requirements
        mkdir -p $HERE/load
        cd $HERE/load
        TINI_VERSION=v0.19.0
        ! [ -f tini ] && wget https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini
        cd $HERE

        export DOCKER_BUILDKIT=1
        docker build \
            --build-arg="CONDA_ENV=${NAME}_env" \
            -t $DOCKER_IMAGE:$VERSION -t $DOCKER_IMAGE:latest .
    ;;
    --push|-p)
        # sudo docker login quay.io
        docker push $DOCKER_IMAGE:$VERSION
        docker push $DOCKER_IMAGE:latest
    ;;
    --sif)
        # apptainer image from the LOCAL docker daemon (no registry round-trip)
        apptainer build $HERE/$NAME.sif docker-daemon://$DOCKER_IMAGE:$VERSION
    ;;
    --check|-c)
        # re-run the image's own pin assertions against the built image
        docker run --rm $DOCKER_IMAGE:$VERSION python -c "\
import numpy, torch, networkx; \
print('numpy', numpy.__version__, '| torch', torch.__version__, \
      '| cuda', torch.version.cuda, '| nx', networkx.__version__)"
    ;;
    --run|-r)
        docker run -it --rm \
            --mount type=bind,source="$HERE",target="/ws" \
            --workdir="/ws" \
            -u $(id -u):$(id -g) \
            $DOCKER_IMAGE:$VERSION \
            /bin/bash
    ;;
    *)
        echo "usage: dev.sh [--build|--push|--sif|--check|--run]"
        echo "bad option: $1"
    ;;
esac
