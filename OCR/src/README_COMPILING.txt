The source is given as a Qt project, though there are no dependencies to any of the Qt libraries.
However, whichever environment is used, the user must ensure that the equivalent to the following Qt
compiler directives are handled : 
QMAKE_CXXFLAGS += -std=c++0xThe project requires the use of c++0x (for the unordered_set and the unordered_map classes).
If Qt is installed, the simplest to to get into the folder, and enter the two commands
> qmake SuperGeneTrees.pro "CONFIG+=RELEASE"
> make

Otherwise, the provided Makefile might be enough to simply run 
> make
And otherwise, the simplest way is to include all .h and .cpp in an Eclipse or QCreator project (or even VC++), and build.