HERE=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

NAME=fabfos
DOCKER_IMAGE=quay.io/hallamlab/$NAME
VER=$(cat $HERE/src/$NAME/version.txt)
# CONDA=conda
CONDA=mamba # https://mamba.readthedocs.io/en/latest/mamba-installation.html#mamba-install
echo image: $DOCKER_IMAGE:$VER
echo ""

# this file contains a list of commands useful for dev,
# providing automation for some build tasks
#
# example workflow 1, pip:
# dev.sh --idev # create a local conda dev env
# # add pypi api token as file to ./secrets [https://pypi.org/help/#apitoken]
# # make some changes to source
# # bump up ./src/fabfos/version.txt
# dev.sh -bp # build the pip package
# dev.sh -up # test upload to testpypi
# dev.sh -upload-pypi # release to pypi index for pip install
#
# example workflow 2, conda:
# dev.sh --idev # create a local conda dev env
# dev.sh -bp # build the pip package
# dev.sh -bc # build conda package from pip package
# dev.sh -uc # publish to conda index
#
# example workflow 3, containerization:
# dev.sh --idev # create a local conda dev env
# dev.sh -bd # build docker image
# dev.sh -ud # publish to quay.io
# dev.sh -bs # build singularity image from local docker image

case $1 in
    ###################################################
    # environments

    --idev) # with dev tools for packaging
        cd $HERE/envs
        $CONDA env create --no-default-packages -n $NAME -f ./base.yml
        $CONDA env update -n $NAME -f ./dev.yml
    ;;
    --ibase) # base only
        cd $HERE/envs
        $CONDA env create --no-default-packages -n $NAME -f ./base.yml
    ;;

    ###################################################
    # build

    -bp) # pip
        # build pip package
        rm -r build
        rm -r dist
        python -m build
    ;;
    -bpi) # pip - test install
        python setup.py install
    ;;
    -bpx) # pip - remove package
        pip uninstall -y $NAME
    ;;
    -bc) # conda
        rm -r $HERE/conda_build
        python ./conda_recipe/compile_recipe.py
        $HERE/conda_recipe/call_build.sh
    ;;
    -bd) # docker
        docker build -t $DOCKER_IMAGE:$VER .
    ;;
    -bs) # singularity image *from docker*
        singularity build $NAME.sif docker-daemon://$DOCKER_IMAGE:$VER
    ;;

    ###################################################
    # upload

    -up) # pip (testpypi)
        PYPI=testpypi
        TOKEN=$(cat secrets/${PYPI}) # https://pypi.org/help/#apitoken
        python -m twine upload --repository $PYPI dist/*.whl -u __token__ -p $TOKEN
    ;;
    -upload-pypi) # pip (pypi)
        echo "not all dependencies are available on pypi, so this is not a good idea..."
        # PYPI=pypi
        # TOKEN=$(cat secrets/${PYPI}) # https://pypi.org/help/#apitoken
        # python -m twine upload --repository $PYPI dist/*.whl -u __token__ -p $TOKEN
    ;;
    -uc) # conda (personal channel)
        # run `anaconda login` first
        find ./conda_build -name *.tar.bz2 | xargs -I % anaconda upload %
    ;;
    -ud) # docker
        # login and push image to quay.io
        # sudo docker login quay.io
	    docker push $DOCKER_IMAGE:$VER
        echo "!!!"
        echo "remember to update the \"latest\" tag"
        echo "https://$DOCKER_IMAGE?tab=tags"
    ;;
    
    ###################################################
    # run

    -r)
        shift
        export PYTHONPATH=$HERE/src:$PYTHONPATH
        python -m $NAME $@
    ;;
    -rd) # docker
            # -e XDG_CACHE_HOME="/ws"\
        shift
        docker run -it --rm \
            -u $(id -u):$(id -g) \
            --mount type=bind,source="$HERE",target="/ws"\
            --workdir="/ws" \
            $DOCKER_IMAGE:$VER /bin/bash
    ;;
    -rs) # singularity
            # -e XDG_CACHE_HOME="/ws"\
        shift
        singularity exec \
            --bind ./:/ws \
            --workdir /ws \
            $HERE/$NAME.sif fabfos /bin/bash
    ;;
    -rt) # single manual test
            # --size 20 \

        export PYTHONPATH=$HERE/src:$PYTHONPATH
        cd scratch
        python -m $NAME \
            --overwrite \
            --threads 12 \
            --output ./no_miffed_test \
            --assembler megahit \
            -i --reads ./beaver_cecum_2ndhits/EKL/Raw_Data/EKL_Cecum_ligninases_pool_secondary_hits_ss01.fastq \
            -b ./ecoli_k12_mg1655.fasta \
            --ends ./beaver_cecum_2ndhits/endseqs.fasta \
            --ends-name-regex "\\w+_\\d+" \
            --ends-fw-flag "FW" \
            --vector ./pcc1.fasta

            # --pool-size 20

            # -i --reads ./beaver_cecum_2ndhits/EKL/Raw_Data/EKL_Cecum_ligninases_pool_secondary_hits_ss10.fastq \
            # --nanopore_reads beaver_cecum_2ndhits/EKL/Raw_Data/EKL_Cecum_ligninases_pool_secondary_hits_ss01.fastq \
    ;;
    *)
        echo "bad option"
        echo $1
    ;;
esac
