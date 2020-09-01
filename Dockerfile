FROM kbase/sdkbase2:python
MAINTAINER nkkchem@gmail.com
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

ENV     NWCHEM_TOP="/codes/nwchem-7.0.0"  
ENV     NWCHEM_DATA=${NWCHEM_TOP}/data 

ENV     USE_MPI=y \
        USE_MPIF=y \
        USE_MPIF4=y \
        NWCHEM_TARGET=LINUX64 \
        BLASOPT="-lopenblas -lpthread -lrt" \
	LAPACK_LIB="-lopenblas -lpthread -lrt" \
        BLAS_SIZE=4 \
        USE_64TO32=y \
        NWCHEM_MODULES="smallqm" \
        NWCHEM_EXECUTABLE=${NWCHEM_TOP}/bin/LINUX64/nwchem \
        NWCHEM_BASIS_LIBRARY=${NWCHEM_DATA}/libraries/ \
        NWCHEM_NWPW_LIBRARY=${NWCHEM_DATA}/libraryps/  \
        FFIELD=amber  \
        AMBER_1=${NWCHEM_DATA}/amber_s/  \
        AMBER_2=${NWCHEM_DATA}/amber_q/  \
        AMBER_3=${NWCHEM_DATA}/amber_x/  \
        AMBER_4=${NWCHEM_DATA}/amber_u/  \
        SPCE=${NWCHEM_DATA}/solvents/spce.rst  \
        CHARMM_S=${NWCHEM_DATA}/charmm_s/  \
        CHARMM_X=${NWCHEM_DATA}/charmm_x/ 

#Build NWChem

RUN     apt-get update \
        && apt-get -y upgrade \
        && apt-get -y install wget bzip2 ssh python-dev gfortran libopenblas-dev libopenmpi-dev \
        openmpi-bin tcsh make openbabel \
        &&  apt-get clean \
        && mkdir codes \
        && cd codes \
        && wget https://github.com/nwchemgit/nwchem/releases/download/v7.0.0-release/nwchem-7.0.0-release.revision-2c9a1c7c-src.2020-02-26.tar.bz2 \
	&& tar -vxjf nwchem-7.0.0-release.revision-2c9a1c7c-src.2020-02-26.tar.bz2 \
        && rm nwchem-7.0.0-release.revision-2c9a1c7c-src.2020-02-26.tar.bz2 \
        && cd ${NWCHEM_TOP}/src \
        && make nwchem_config && make 64_to_32  && make -j3 && \
        mkdir ${NWCHEM_DATA} && \
        mv ${NWCHEM_TOP}/src/basis/libraries ${NWCHEM_DATA} && \
        mv ${NWCHEM_TOP}/src/nwpw/libraryps ${NWCHEM_DATA} && \
        mv ${NWCHEM_TOP}/src/data/amber_s ${NWCHEM_DATA} && \
        mv ${NWCHEM_TOP}/src/data/amber_q ${NWCHEM_DATA} && \
        mv ${NWCHEM_TOP}/src/data/amber_x ${NWCHEM_DATA} && \
        mv ${NWCHEM_TOP}/src/data/amber_u ${NWCHEM_DATA} && \
        mv ${NWCHEM_TOP}/src/data/solvents ${NWCHEM_DATA} && \
        mv ${NWCHEM_TOP}/src/data/charmm_s ${NWCHEM_DATA} && \
        mv ${NWCHEM_TOP}/src/data/charmm_x ${NWCHEM_DATA} && \
        rm -rf $NWCHEM_TOP/src && \
        rm -rf $NWCHEM_TOP/lib 

RUN     apt-get -y remove  ssh tcsh  gfortran  python-dev libopenmpi-dev && apt-get clean
RUN     conda update -n base -c defaults conda
RUN	pip install snakemake
#RUN     conda install -c conda-forge -c bioconda snakemake
RUN     conda install -c rdkit rdkit
RUN     conda install -c openbabel openbabel

ENV     NWCHEM_SIM_DIR="/simulation"
ENV     NWCHEM_BIN=${NWCHEM_TOP}/bin/LINUX64
ENV     NWCHEM_TEMPLATES_DIR=${NWCHEM_DATA}/templates
ENV     PATH="${NWCHEM_BIN}:$PATH"

COPY ./nwchem-scripts/inchi_to_submission.py ${NWCHEM_BIN}/
COPY ./nwchem-scripts/extract_properties_mulliken_charges_mol2.py ${NWCHEM_BIN}/
COPY ./nwchem-scripts/compound_parsing.py ${NWCHEM_BIN}/
COPY ./nwchem-scripts/export.py ${NWCHEM_BIN}/
COPY ./nwchem-scripts/CompoundSetUtil.py ${NWCHEM_BIN}/
COPY ./nwchem-scripts/extract_properties.py ${NWCHEM_BIN}/
COPY ./nwchem-scripts/inchiTo3D.py ${NWCHEM_BIN}/
COPY ./nwchem-scripts/dft-cme.py ${NWCHEM_BIN}/
COPY ./nwchem-scripts/create_nw_files.py ${NWCHEM_BIN}
COPY ./snakemake-scripts/final_pipeline.snakemake ${NWCHEM_BIN}
COPY ./snakemake-scripts/cme.py ${NWCHEM_BIN}
COPY ./dft_templates/template_dft.nw ${NWCHEM_BIN}
COPY ./dft_templates/template_dft.sbatch ${NWCHEM_BIN}
COPY ./dft_templates/template_odft.nw ${NWCHEM_BIN}
RUN   mkdir ${NWCHEM_SIM_DIR}

# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

#ENTRYPOINT [ "./scripts/entrypoint.sh" ]
#CMD [ ]
