#ifndef SHADER_H
#define SHADER_H
// =====================================================================
// This file is part of SWE_CUDA (see file SWE_Block.cu for details).
// 
// Copyright (C) 2010,2011 Tobias Schnabel
// 
// SWE_CUDA is free software: you can redristribute it and/or modify
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
#include <SDL.h>
#include <SDL_opengl.h>
#include <iostream>
#include <fstream>

class Shader {
public:
	// Constructor and Destructor
	// Load vertex shader file and fragment shader file 
	Shader(char const * vertexShaderFile,  char const * fragmentShaderFile) ;
	~Shader();

	// Check if shaders are supported
	bool shadersSupported();
	bool shadersLoaded();

	// Shader control
	void enableShader();
	void disableShader();

private:	
	// State flags
	bool shdrSupport;
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

	// Shader extension function pointers
	PFNGLCREATESHADERPROC glCreateShader;
	PFNGLCREATEPROGRAMPROC glCreateProgram;
	PFNGLATTACHSHADERPROC glAttachShader;
	PFNGLCOMPILESHADERPROC glCompileShader;
	PFNGLUSEPROGRAMPROC glUseProgram;
	PFNGLDETACHSHADERPROC glDetachShader;
	PFNGLDELETESHADERPROC glDeleteShader;
	PFNGLLINKPROGRAMPROC glLinkProgram;
	PFNGLSHADERSOURCEPROC glShaderSource;
	PFNGLDELETEPROGRAMPROC glDeleteProgram;
	
	// Shader objects extension pointers
	PFNGLGETOBJECTPARAMETERIVARBPROC glGetObjectParameterivARB;
	PFNGLGETSHADERIVPROC glGetShaderiv;
	PFNGLGETSHADERINFOLOGPROC glGetShaderInfoLog;
	PFNGLGETPROGRAMINFOLOGPROC glGetProgramInfoLog;

};
#endif
