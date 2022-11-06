SWE
===

[An Education-Oriented Code for Parallel Tsunami Simulation.](Documentation/Mainpage.md)

## General Information
This code has been tested under Ubuntu 20.04 or higher. Other Linux distributions will probably work. However, we do not recommend using Windows or MacOS. If you do not have Linux installed on your computer, please use WSL, Virtual Machine, or similar.
Since SWE1D uses CMake, it should work under Windows with MSVC or MacOS.

## Dockerfile
You can use the Dockerfile included in the repository. This will create an image where everything is set up and you can just begin with compiling the code.

Of course, you can also just look at the Dockerfile as a recipe on how to set up your environment.

For a short introduction on how to generate, build and run a Docker container, refer to [Dockerfile](/Tools/README.md).

## Tutorial
### Compilation
As build system configurator we use CMake. To compile the code execute the following commands in this directory:

* Create a build directory: `mkdir build`. You can also choose any other name for your build directory.
* Switch to this directory: `cd build`
* (Optional): Choose the compiler being used (if you want to use a specific MPI compiler/version): `export CXX=mpic++`
* Run CMake: `cmake ..` (this configures a `RelWithDebInfo` build, which is default. For a `Debug` build, run `cmake .. -DCMAKE_BUILD_TYPE=Debug`) and for a `Release` build, `cmake .. -DCMAKE_BUILD_TYPE=Release`. This is especially recommended in production and benchmark runs.
* For developing, consider `cmake .. -DENABLE_DEVELOPER_MODE=ON`. For an overview of all availble options, use `ccmake ..`
* Run Make: `make` (or `make -j` for compiling with multiple cores)
* Run `make help` to see all available targets to build.

### Running a Simulation
* Run the code in serial via `./SWE-Runner`
* With `./SWE-Runner --help`, you can see additional command-line arguments you can pass.

### Adding new source files
You can add new source files by just creating them somewhere within the `Source` folder. CMake automatically detects these files and adds them to the build.

### Testing
Some basic unit tests have been implemented (`make test`). Feel free to add your own test cases inside the `Tests` folder.

### Creating a Doxygen documentation
* Run the following CMake command: `cmake .. -DENABLE_DEVELOPER_MODE=ON -DOPT_ENABLE_DOXYGEN=ON`
* `make doxygen-docs`
* You should now see in your build directory a folder called `html`
* Open `index.html` to see the generated documentation.

## Development Hints & FAQ
### It does not compile and everything seems fine?
Make sure to use `make clean` before you use `make`. Sometimes there are build artifacts from previous build processes that spoil your current compilation process. `make clean` takes care of deleting everything that should not be there and allows the compiler to start from scratch.

Sometimes it is also helpful to delete the `build` folder and create a new one, following the steps from the compilation section above.

### How can I see all the compiler flags the generated Makefile is using?
Instead of using `make`, run `VERBOSE=1 make`. You can also run `make -n` to invoke a dry run where you see what the Makefile would do in case of compilation.

### How can I see the test output?
Instead of using `make test`, run `ctest --verbose`.