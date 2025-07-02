
# Config for OLCF Frontier

if [[ $HOST == *".aurora.alcf.anl.gov" ]]
then
  # We seem to do fine without this
  # HOST_ARCH=
  # DEVICE_ARCH=

  MPI_EXE=mpiexec
  NPROC=32

  module load cmake hdf5

  EXTRA_FLAGS="-DPARTHENON_DISABLE_HDF5_COMPRESSION=ON"

  if [[ $ARGS == *"sycl"* ]]; then
    # SYCL compile for Intel GPUs

    #EXTRA_FLAGS="-DKHARMA_SPLIT_IMPLICIT_SOLVE=ON $EXTRA_FLAGS"

    # CXX_NATIVE=CC
    # C_NATIVE=cc

    # Runtime: WTF Intel, this is so complicated
    export CPU_BIND_SCHEME="--cpu-bind=list:1-8:9-16:17-24:25-32:33-40:41-48:53-60:61-68:69-76:77-84:85-92:93-100"
    MPI_NUM_PROCS=${MPI_NUM_PROCS:-12}
    MPI_EXTRA_ARGS="-ppn 12 $CPU_BIND_SCHEME gpu_tile_compact.sh"

    export OMP_PROC_BIND=${OMP_PROC_BIND:-spread}
    export OMP_PLACES=${OMP_PLACES:-threads}

  else
    if [[ $ARGS == *"gcc"* ]]; then
      # GCC already loaded, just specify
      CXX_NATIVE=g++
      C_NATIVE=gcc
    fi

    # Runtime
    MPI_NUM_PROCS=1

    export OMP_PROC_BIND=${OMP_PROC_BIND:-spread}
    export OMP_PLACES=${OMP_PLACES:-threads}
  fi
fi

if [[ $HOST == *".crux.alcf.anl.gov" ]]
then
  HOST_ARCH=ZEN3

  MPI_EXE=mpiexec

  # """load""" CMake. C'mon guys
  export PATH=$HOME/libs/bin:$PATH

  if [[ $ARGS == *"gcc"* ]]; then
    module swap PrgEnv-cray PrgEnv-gnu
  fi
  module load cray-hdf5-parallel
  
  EXTRA_FLAGS="-DPARTHENON_DISABLE_HDF5_COMPRESSION=ON"
  
  C_NATIVE=cc
  CXX_NATIVE=CC

  # CPU Compile
  # TODO -c etc etc
  MPI_NUM_PROCS=1
fi