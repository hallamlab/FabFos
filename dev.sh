NAME=fabfos
DOCKER_IMAGE=quay.io/hallam_lab/$NAME
VER=1.9
echo image: $DOCKER_IMAGE:$VER
echo ""

HERE=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

case $1 in
    --build|-b)
        # change the url in python if not txyliu
        # build the docker container locally *with the cog db* (see above)
        docker build -t $DOCKER_IMAGE:$VER .
    ;;
    --push|-p)
        # login and push image to quay.io
        # sudo docker login quay.io
	    docker push $DOCKER_IMAGE:$VER
    ;;
    --sif)
        # test build singularity
        singularity build $NAME.sif docker-daemon://$DOCKER_IMAGE:$VER
    ;;
    --run|-r)
        # test run docker image
            # --mount type=bind,source="$HERE/scratch",target="/ws" \
            # --mount type=bind,source="$HERE/scratch/res",target="/ref"\
            # -e XDG_CACHE_HOME="/ws"\
            # --workdir="/ws" \
            # -u $(id -u):$(id -g) \
        shift
        docker run -it --rm $DOCKER_IMAGE fabfos $@

    ;;
    
    -t)
        #
        # scratch space for testing stuff
        #
        cd scratch
        docker run -it --rm \
            --mount type=bind,source="./",target="/ws" \
            --workdir="/ws" \
            -u $(id -u):$(id -g) \
            $DOCKER_IMAGE fabfos \
            --threads 14 --fabfos_path ./ --force --assembler megahit \
            -m miffed.csv --reads ./reads -i -b ecoli_k12_mg1655.fasta
    ;;
    *)
        echo "bad option"
        echo $1
    ;;
esac
