#!/bin/bash

# Astyle required http://astyle.sourceforge.net/
astyle --style=allman --recursive --suffix=none 'include/*.h'
astyle --style=allman --recursive --suffix=none 'include/*.h.in'
astyle --style=allman --recursive --suffix=none 'src/*.c'
astyle --style=allman --recursive --suffix=none 'test/*.c'
astyle --style=allman --recursive --suffix=none 'examples/*.c'
