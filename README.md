# dMC — differentiable Muskingum-Cunge in C++

dMC is a C++ implementation of a differentiable Muskingum-Cunge channel-routing scheme — suitable for hydrological routing and coupling with models like SUMMA.

## Features

- C++ implementation with clean modular structure (include / src / external)  
- Build via CMake  
- NetCDF support 
- Example executables & test suite included  
- GPL-3.0 license  

## Requirements / Dependencies

- CMake (>= 3.15 recommended)  
- A C++ compiler with C++17 support  
- NetCDF library 

## Build Instructions

```bash
git clone https://github.com/DarriEy/dMC.git
cd dMC
mkdir build
cd build
cmake .. -DDMC_ENABLE_NETCDF=ON    # add other flags as needed
make -j4
```

## Running an Example

```bash
./dmc_route_run \
  -n /path/to/topology.nc \
  -f /path/to/run_1_timestep.nc \
  -p summa \
  -o discharge.csv \
  -j jacobian.csv
```

Or if you want to run with config.txt

```bash
./dmc_route_run \
  -c /path/to/config.txt \
```

## Project Structure

```
dMC/
  ├── external/       # external dependencies (e.g. BMI, NetCDF wrappers, etc.)
  ├── include/        # public headers
  ├── src/            # implementation code
  ├── tests/          # unit or regression tests
  ├── examples/       # runnable example configurations or data
  ├── CMakeLists.txt  # main build configuration
  ├── DESIGN.md       # design / architecture notes
  ├── LICENSE         # GPL-3.0 license
  └── README.md       # this file
```

## How to Contribute

1. Fork the repository.  
2. Create a branch for your feature / fix.  
3. Run tests and ensure everything compiles.  
4. Submit a pull request with a clear explanation of changes.  

## License

This project is licensed under the GPL-3.0 license — see the `LICENSE` file for details.

