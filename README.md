# dMC
differentiable Muskingum Cunge in C++

## Build instructions

git clone https://github.com/DarriEy/dMC.git
cd dMC
mkdir build
cd build
cmake .. -DDMC_ENABLE_NETCDF=ON
make -j4


test run ./dmc_route_run

## Usage instructions (for SUMMA with jacobian)

./dmc_route_run \
  -n /path/to/topology.nc \
  -f /path/to/run_1_timestep.nc \
  -p summa \
  -o discharge.csv \
  -j jacobian.csv