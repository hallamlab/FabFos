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

    ###################################################
    # test

    -t) # single manual test
            # --size 20 \

        # --overwrite \
            # -a spades_meta \
            # --endf ../data/beaver/endseq_CEC_FW.fa \
            # --endr ../data/beaver/endseq_CEC_RE.fa \
            # --endf ../data/beaver/endseq_COL_FW.fa \
            # --endr ../data/beaver/endseq_COL_RE.fa \
            # -s ./inputs/ss01.fastq \
            # -i ../data/beaver/Beaver_colon/2nd_hits/EKL/Raw_Data/EKL_Colon_ligninases_pool_secondary_hits.fastq ../data/beaver/Beaver_colon/2nd_hits/EOL/Raw_Data/EOL_Colon_ligninases_pool_secondary_hits.fastq \
            # -i ../data/beaver/Beaver_colon/2nd_hits/EKL/Raw_Data/EKL_Colon_ligninases_pool_secondary_hits.fastq \
            # -a megahit \
            # --overwrite \
            # -a megahit_sensitive megahit_meta spades_isolate spades_meta \
        # /home/tony/workspace/grad/tools/fabfos/scratch/test01_assemblies/spades_isolate/contigs.fasta
        # /home/tony/workspace/grad/tools/fabfos/scratch/test01_assemblies/spades_isolate/contigs.fasta
        # /home/tony/workspace/grad/tools/FabFos/scratch/test01_assemblies/spades_isolate/contigs.fasta
            # -b ./ecoli_k12_mg1655.fasta \
            # -a \
            # /home/tony/workspace/grad/tools/FabFos/data/beaver/112.bvr_compound_assembly.2023-11-22-23-18/cec_1/temp_assembly/megahit_sensitive/final.contigs.fa \
            # /home/tony/workspace/grad/tools/FabFos/data/beaver/112.bvr_compound_assembly.2023-11-22-23-18/cec_1/temp_assembly/megahit_meta/final.contigs.fa \
            # /home/tony/workspace/grad/tools/FabFos/data/beaver/112.bvr_compound_assembly.2023-11-22-23-18/cec_1/temp_assembly/spades_isolate/contigs.fasta \
            # /home/tony/workspace/grad/tools/FabFos/data/beaver/112.bvr_compound_assembly.2023-11-22-23-18/cec_1/temp_assembly/spades_meta/contigs.fasta \
            # /home/tony/workspace/grad/tools/FabFos/data/beaver/112.bvr_compound_assembly.2023-11-22-23-18/cec_2/temp_assembly/megahit_sensitive/final.contigs.fa \
            # /home/tony/workspace/grad/tools/FabFos/data/beaver/112.bvr_compound_assembly.2023-11-22-23-18/cec_2/temp_assembly/megahit_meta/final.contigs.fa \
            # --endf ../data/beaver/endseq_CEC_FW.fa \
            # --endr ../data/beaver/endseq_CEC_RE.fa \
            # -a megahit_sensitive \
            # ../data/beaver/112.bvr_compound_assembly.2023-11-22-23-18/cec_2/temp_assembly/spades_isolate/contigs.fasta \
            # ../data/beaver/112.bvr_compound_assembly.2023-11-22-23-18/cec_2/temp_assembly/spades_meta/contigs.fasta \
            # ../data/beaver/112.bvr_compound_assembly.2023-11-22-23-18/cec_2/temp_assembly/spades_meta/contigs.fasta \
            # --end_regex "\w+_\d+" \
        export PYTHONPATH=$HERE/src:$PYTHONPATH
        cd scratch
        python -m $NAME run -t 12 \
            -i ./inputs/ss10.fastq.gz \
            -d ../../../resources/mp3db \
            --vector ./inputs/pcc1.fasta \
            -o ./test01
    ;;

    -tm)
    cd scratch
    metapathways run -t 12 \
        -i ./mpwi/contigs.fna \
        -r ./mpwi/reads \
        -d ../../../resources/mp3db \
        --annotation_dbs swissprot metacyc\
        --COMPUTE_TPM skip \
        -o ./mpw_test
    ;;

    -td) # test docker
        shift
        docker run -it --rm \
            -u $(id -u):$(id -g) \
            --mount type=bind,source="$HERE",target="/ws"\
            --mount type=bind,source="$HERE/src/fabfos",target="/app/fabfos"\
            --mount type=bind,source="/home/tony/workspace/grad/analysis/beaver_microbiome/data/fabfos",target="/home/tony/workspace/grad/analysis/beaver_microbiome/data/fabfos"\
            --workdir="/ws" \
            $DOCKER_IMAGE:$VER fabfos run -t 12 \
                -a megahit \
                --overwrite \
                -i /ws/scratch/inputs/ss01.fastq.gz \
                -b /ws/scratch/ecoli_k12_mg1655.fasta \
                --endf /ws/data/beaver/endseq_COL_FW.fa \
                --endr /ws/data/beaver/endseq_COL_RE.fa \
                --end_regex "\\w+_\\d+" \
                -o /ws/scratch/test_docker01
    ;;

    -ts) # test singularity
        shift
        singularity run -B $HERE:/ws,"/home/tony/workspace/grad/analysis/beaver_microbiome/data/fabfos":"/home/tony/workspace/grad/analysis/beaver_microbiome/data/fabfos" \
        $HERE/fabfos.sif fabfos run -t 12 \
            -a megahit \
            --overwrite \
            -i /ws/scratch/inputs/ss10.fastq.gz \
            -b /ws/scratch/ecoli_k12_mg1655.fasta \
            --endf /ws/data/beaver/endseq_COL_FW.fa \
            --endr /ws/data/beaver/endseq_COL_RE.fa \
            --end_regex "\\w+_\\d+" \
            -o /ws/scratch/test_singularity01
    ;;

    ###################################################

    *)
        echo "bad option"
        echo $1
    ;;
esac
