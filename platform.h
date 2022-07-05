#if !defined(PLATFORM_LINUX) && !defined(PLATFORM_WINDOWS)
#define PLATFORM_LINUX
#endif

#if defined(PLATFORM_LINUX)
#include "platformlinux.h"
#include "SDL2/SDL.h"
#else // PLATFORM_WINDOWS
#include "platformwin.h"
#include "SDL.h"
#endif

#include <stdio.h>
