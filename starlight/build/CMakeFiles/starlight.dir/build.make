# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.26

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/xihe/starlight

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/xihe/starlight/build

# Include any dependencies generated for this target.
include CMakeFiles/starlight.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/starlight.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/starlight.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/starlight.dir/flags.make

CMakeFiles/starlight.dir/src/main.cpp.o: CMakeFiles/starlight.dir/flags.make
CMakeFiles/starlight.dir/src/main.cpp.o: /home/xihe/starlight/src/main.cpp
CMakeFiles/starlight.dir/src/main.cpp.o: CMakeFiles/starlight.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/xihe/starlight/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/starlight.dir/src/main.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/starlight.dir/src/main.cpp.o -MF CMakeFiles/starlight.dir/src/main.cpp.o.d -o CMakeFiles/starlight.dir/src/main.cpp.o -c /home/xihe/starlight/src/main.cpp

CMakeFiles/starlight.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/starlight.dir/src/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/xihe/starlight/src/main.cpp > CMakeFiles/starlight.dir/src/main.cpp.i

CMakeFiles/starlight.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/starlight.dir/src/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/xihe/starlight/src/main.cpp -o CMakeFiles/starlight.dir/src/main.cpp.s

# Object files for target starlight
starlight_OBJECTS = \
"CMakeFiles/starlight.dir/src/main.cpp.o"

# External object files for target starlight
starlight_EXTERNAL_OBJECTS =

starlight: CMakeFiles/starlight.dir/src/main.cpp.o
starlight: CMakeFiles/starlight.dir/build.make
starlight: libStarlib.a
starlight: CMakeFiles/starlight.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/xihe/starlight/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable starlight"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/starlight.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/starlight.dir/build: starlight
.PHONY : CMakeFiles/starlight.dir/build

CMakeFiles/starlight.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/starlight.dir/cmake_clean.cmake
.PHONY : CMakeFiles/starlight.dir/clean

CMakeFiles/starlight.dir/depend:
	cd /home/xihe/starlight/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/xihe/starlight /home/xihe/starlight /home/xihe/starlight/build /home/xihe/starlight/build /home/xihe/starlight/build/CMakeFiles/starlight.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/starlight.dir/depend

