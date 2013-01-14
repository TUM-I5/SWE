#ifndef SHADER_H
#define SHADER_H
// =====================================================================
// This file is part of SWE_CUDA (see file SWE_Block.cu for details).
// 
// Copyright (C) 2010,2011 Tobias Schnabel
// Copyright (C) 2012      Sebastian Rettenberger
// 
// SWE_CUDA is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// SWE_CUDA is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with SWE_CUDA.  If not, see <http://www.gnu.org/licenses/>.
// =====================================================================
#include <SDL/SDL.h>
#include <SDL/SDL_opengl.h>
#include <iostream>
#include <fstream>

class Shader {
public:
	// Constructor and Destructor
	// Load vertex shader file and fragment shader file 
	Shader(char const * vertexShaderFile,  char const * fragmentShaderFile) ;
	~Shader();

	// Check if shaders are loaded
	bool shadersLoaded();

	// Shader control
	void enableShader();
	void disableShader();

	/**
	 * @return Location of the uniform variable
	 */
	GLint getUniformLocation(const char* name)
	{
		if (!shdrLoaded)
			return -1;
		return glGetUniformLocation(program, name);
	}

	/**
	 * Set a uniform variable in the shader
	 */
	void setUniform(GLint location, GLfloat value)
	{
		if (location < 0)
			return;
		glUniform1f(location, value);
	}

private:	
	// State flags
	bool shdrLoaded;

	// Helper functions
	bool readShaderFile(char const * filename, GLchar * & shaderSource, GLint & length);
	bool isExtensionSupported(const char* szTargetExtension );
	bool isShaderCompiled(GLuint shader, char const * prefix);
	bool isProgramLinked(GLuint program, char const * prefix);

	// Vertex-Shader
	GLuint   vertexShader;
	GLchar * vertexShaderSource;
	GLint    vertexShaderLength;
	
	// Fragment-Shader
	GLuint   fragmentShader;
	GLchar * fragmentShaderSource;
	GLint    fragmentShaderLength;

	// Shaders id
	GLuint program;

	/** Are shaders supported */
	static bool shdrSupport;

	// Shader extension function pointers
	static PFNGLCREATESHADERPROC glCreateShader;
	static PFNGLCREATEPROGRAMPROC glCreateProgram;
	static PFNGLATTACHSHADERPROC glAttachShader;
	static PFNGLCOMPILESHADERPROC glCompileShader;
	static PFNGLUSEPROGRAMPROC glUseProgram;
	static PFNGLDETACHSHADERPROC glDetachShader;
	static PFNGLDELETESHADERPROC glDeleteShader;
	static PFNGLLINKPROGRAMPROC glLinkProgram;
	static PFNGLSHADERSOURCEPROC glShaderSource;
	static PFNGLDELETEPROGRAMPROC glDeleteProgram;
	static PFNGLGETUNIFORMLOCATIONPROC glGetUniformLocation;
	static PFNGLUNIFORM1FPROC glUniform1f;
	
	// Shader objects extension pointers
	static PFNGLGETOBJECTPARAMETERIVARBPROC glGetObjectParameterivARB;
	static PFNGLGETSHADERIVPROC glGetShaderiv;
	static PFNGLGETSHADERINFOLOGPROC glGetShaderInfoLog;
	static PFNGLGETPROGRAMINFOLOGPROC glGetProgramInfoLog;
};
#endif
