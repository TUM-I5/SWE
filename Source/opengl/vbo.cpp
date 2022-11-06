/**
 * @file
 * This file is part of SWE.
 *
 * @author Sebastian Rettenberger (rettenbs AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
 *
 * @section LICENSE
 *
 * SWE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SWE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SWE.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "vbo.h"
#include "visualization.h"

void VBO::init()
{
	if (glGenBuffers == 0L) {
		// Load vbo extension

		// Check OpenGL extension(s)
		if (!Visualization::isExtensionSupported("GL_ARB_vertex_buffer_object")) {
			tools::Logger::logger.printString("Vertex Buffer Objects Extension not supported! Exit..\n");
			SDL_Quit();
			exit(1);
		}

		// Load Vertex Buffer Extension
		glGenBuffers = (PFNGLGENBUFFERSARBPROC) SDL_GL_GetProcAddress("glGenBuffersARB");
		glBindBuffer = (PFNGLBINDBUFFERARBPROC) SDL_GL_GetProcAddress("glBindBufferARB");
		glBufferData = (PFNGLBUFFERDATAARBPROC) SDL_GL_GetProcAddress("glBufferDataARB");
		glDeleteBuffers = (PFNGLDELETEBUFFERSARBPROC) SDL_GL_GetProcAddress("glDeleteBuffersARB");
	}

	glGenBuffers(1, &name);
}

PFNGLGENBUFFERSARBPROC VBO::glGenBuffers = 0L;
PFNGLBINDBUFFERARBPROC VBO::glBindBuffer = 0L;
PFNGLBUFFERDATAARBPROC VBO::glBufferData = 0L;
PFNGLDELETEBUFFERSARBPROC VBO::glDeleteBuffers = 0L;
