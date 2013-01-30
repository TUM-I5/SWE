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
 *
 *
 * @section DESCRIPTION
 *
 * Handles a VertexBufferObject.
 */
#ifndef VBO_H
#define VBO_H

#include "tools/Logger.hh"

#include <SDL/SDL_opengl.h>

class VBO
{
private:
	/** OpenGL name of the object */
	GLuint name;

public:
	VBO()
		: name(0)
	{}

	/**
	 * Initializes the object
	 */
	void init();

	/**
	 * @return The OpenGL name of the buffer
	 */
	GLuint getName()
	{
		return name;
	}

	void setBufferData(GLsizei size, const void* data,
			GLenum target = GL_ARRAY_BUFFER,
			GLenum usage = GL_STATIC_DRAW)
	{
		glBindBuffer(target, name);
		glBufferData(target, size, data, usage);
		glBindBuffer(target, 0);
	}

	void bindBuffer(GLenum target = GL_ARRAY_BUFFER)
	{
		glBindBuffer(target, name);
	}

	/**
	 * Frees all associated memory
	 */
	void finialize()
	{
		if (name) {
			glDeleteBuffers(1, &name);
			name = 0;
		}
	}

private:
	// VBO Extension Function Pointers
	static PFNGLGENBUFFERSARBPROC glGenBuffers;					// VBO Name Generation Procedure
	static PFNGLBINDBUFFERARBPROC glBindBuffer;					// VBO Bind Procedure
	static PFNGLBUFFERDATAARBPROC glBufferData;					// VBO Data Loading Procedure
	static PFNGLDELETEBUFFERSARBPROC glDeleteBuffers;			// VBO Deletion Procedure
};

#endif // VBO_H
