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

#include "shader.h"

#include "tools/Logger.hh"

/**
    Constructor.
	Check whether shaders are supported. If yes, load vertex and fragment 
	shader from textfile into memory and compile

	@param  vertexShaderFile	name of the text file containing the vertex 
								shader code
	@param	fragmentShaderFile  name of the text file containing the fragment
								shader code
	
*/	
Shader::Shader(char const * vertexShaderFile, char const * fragmentShaderFile)
{
	if (!shdrSupport) {
		// Shaders are either not supported or we did not try yet
		// This if-statement makes sure, we only load the extensions once

		shdrSupport = isExtensionSupported("GL_ARB_vertex_shader")
						&& isExtensionSupported("GL_ARB_shader_objects")
						&& isExtensionSupported("GL_ARB_fragment_shader");

		if (shdrSupport) {
			tools::Logger::logger.printString("Shaders supported!");

			// Load extensions
			glCreateShader = (PFNGLCREATESHADERPROC) SDL_GL_GetProcAddress("glCreateShader");
			glCreateProgram = (PFNGLCREATEPROGRAMPROC) SDL_GL_GetProcAddress("glCreateProgram");
			glAttachShader = (PFNGLATTACHSHADERPROC)SDL_GL_GetProcAddress("glAttachShader");
			glCompileShader = (PFNGLCOMPILESHADERPROC)SDL_GL_GetProcAddress("glCompileShader");
			glUseProgram = (PFNGLUSEPROGRAMPROC)SDL_GL_GetProcAddress("glUseProgram");
			glDetachShader = (PFNGLDETACHSHADERPROC)SDL_GL_GetProcAddress("glDetachShader");
			glDeleteShader = (PFNGLDELETESHADERPROC) SDL_GL_GetProcAddress("glDeleteShader");
			glLinkProgram = (PFNGLLINKPROGRAMPROC)SDL_GL_GetProcAddress("glLinkProgram");
			glShaderSource = (PFNGLSHADERSOURCEPROC)SDL_GL_GetProcAddress("glShaderSource");
			glDeleteProgram = (PFNGLDELETEPROGRAMPROC)SDL_GL_GetProcAddress("glDeleteProgram");
			glGetUniformLocation = (PFNGLGETUNIFORMLOCATIONPROC) SDL_GL_GetProcAddress("glGetUniformLocation");
			glUniform1f = (PFNGLUNIFORM1FPROC) SDL_GL_GetProcAddress("glUniform1f");
			glGetObjectParameterivARB = (PFNGLGETOBJECTPARAMETERIVARBPROC) SDL_GL_GetProcAddress("glGetObjectParameterivARB");
			glGetShaderiv = (PFNGLGETSHADERIVPROC) SDL_GL_GetProcAddress("glGetShaderiv");
			glGetShaderInfoLog = (PFNGLGETSHADERINFOLOGPROC) SDL_GL_GetProcAddress("glGetShaderInfoLog");
			glGetProgramInfoLog = (PFNGLGETPROGRAMINFOLOGPROC) SDL_GL_GetProcAddress("glGetProgramInfoLog");
		} else
			tools::Logger::logger.printString("Shaders are NOT supported! Normal rendering mode");
	}

	shdrLoaded = false;

	if (shdrSupport) {
		// Read shader files
		bool readSuccess = false;
		vertexShaderSource = NULL;
		fragmentShaderSource = NULL;
		readSuccess = readShaderFile(vertexShaderFile, vertexShaderSource, vertexShaderLength);
		readSuccess &= readShaderFile(fragmentShaderFile, fragmentShaderSource, fragmentShaderLength);

		// Check if reading was completed
		if (readSuccess) {
			// Create Vertex-Shader 
			vertexShader = glCreateShader(GL_VERTEX_SHADER);
			glShaderSource(vertexShader, 1, const_cast<GLchar const **>(&vertexShaderSource), &vertexShaderLength);
			glCompileShader(vertexShader);

			// Create Fragment-Shader 
			fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
			glShaderSource(fragmentShader, 1, const_cast<GLchar const **>(&fragmentShaderSource), &fragmentShaderLength);
			glCompileShader(fragmentShader);

			if (isShaderCompiled(vertexShader, "Vertex-Shader")
				&& isShaderCompiled(fragmentShader, "Fragment-Shader"))
			{
				// Create Shader program
				program = glCreateProgram();
				glAttachShader(program, vertexShader);
				glAttachShader(program, fragmentShader);
				glLinkProgram(program);
				shdrLoaded = true;
			} else {
				// Errors while compiling shaders
				glDeleteShader(vertexShader);
				delete[] vertexShaderSource;
				glDeleteShader(fragmentShader);
				delete[] fragmentShaderSource;
			}
		} else {
			if (vertexShaderSource != NULL) delete[] vertexShaderSource;
			if (fragmentShaderSource != NULL) delete[] fragmentShaderSource;
		}
	}
}
/**
    Destructor.
	Unload shaders and free resources.
	
*/	
Shader::~Shader() {
	// Unload Shaders
	if (shdrLoaded) {
		glUseProgram(0);
		glDeleteProgram(program);

		glDetachShader(program, vertexShader);
		glDeleteShader(vertexShader);
		delete[] vertexShaderSource;

		glDetachShader(program, fragmentShader);
		glDeleteShader(fragmentShader);
		delete[] fragmentShaderSource;
	}
}
/**
    Replaces OpenGL shaders by our custom shaders

*/	
void Shader::enableShader() {
	if (shdrLoaded) {
		glUseProgram(program);
	}
}
/**
    Restores OpenGL default shaders 

*/
void Shader::disableShader() {
	if (shdrLoaded) {
		glUseProgram(0);
	}
}

/**
    Returns, whether shaders could by loaded successfully 

*/
bool Shader::shadersLoaded() {
	return shdrLoaded;
}

/**
    Returns, whether a special extension is supported by the current 
	graphics card

    @param	szTargetExtention	string describing the extension to look for	

*/	
bool Shader::isExtensionSupported(const char* szTargetExtension )
{
	const unsigned char *pszExtensions = NULL;
	const unsigned char *pszStart;
	unsigned char *pszWhere, *pszTerminator;
	
	// Extension names should not have spaces
	pszWhere = (unsigned char *) strchr( szTargetExtension, ' ' );
	if( pszWhere || *szTargetExtension == '\0' )
		return false;

	// Get Extensions String
	pszExtensions = glGetString( GL_EXTENSIONS );

	// Search The Extensions String For An Exact Copy
	pszStart = pszExtensions;
	for(;;)
	{
		pszWhere = (unsigned char *) strstr( (const char *) pszStart, szTargetExtension );
		if( !pszWhere )
			break;
		pszTerminator = pszWhere + strlen( szTargetExtension );
		if( pszWhere == pszStart || *( pszWhere - 1 ) == ' ' )
			if( *pszTerminator == ' ' || *pszTerminator == '\0' )
				return true;
		pszStart = pszTerminator;
	}
	return false;
}

/**
    Returns, whether a special extension is supported by the current 
	graphics card

    @param	filename		shader file
	@return shaderSource 	file content
	@return length			length of the file
	@returns bool	true, if file has been read successfully

*/	
bool Shader::readShaderFile(char const * filename, GLchar * & shaderSource, GLint & length)
{
	using namespace std;
	string fullFileName = "";

	// On unix: determine directory of exec 
#ifdef __unix__
	// Code taken from: 
	// http://www.gamedev.net/community/forums/topic.asp?topic_id=459511
	std::string path = "";
	pid_t pid = getpid();
	char buf[20] = {0};
	sprintf(buf,"%d",pid);
	std::string _link = "/proc/";
	_link.append( buf );
	_link.append( "/exe");
	char proc[512];
	int ch = readlink(_link.c_str(),proc,512);
	if (ch != -1) {
		proc[ch] = 0;
		path = proc;
		std::string::size_type t = path.find_last_of("/");
		path = path.substr(0,t);
	}

	fullFileName = path + string("/");
#endif

	fullFileName.append(filename);
	ifstream file(fullFileName.c_str(), ios::in);

	length = 0;

	if (!file.good())
	{
		return false;
	}

	file.seekg(0, ios::end);
	length = file.tellg();
	file.seekg(ios::beg);

	if (length == 0 )
	{
		return false;
	}

	shaderSource = new GLchar[length];

	if (shaderSource == 0)
	{
		return false;
	}

	memset(shaderSource, 0, length);
	file.read(shaderSource, length);

	file.close();

	return true;
}

/**
    Returns, whether a shader has been compiled with success

    @param	shader		shader id
	@param prefix 		custom error message

*/	
bool Shader::isShaderCompiled(GLuint shader, char const * prefix)
{
	using namespace std;
	GLint compiled(0);

	glGetObjectParameterivARB(shader, GL_COMPILE_STATUS, &compiled);

	if (!compiled)
	{
		GLint maxLength(0);

		glGetShaderiv(shader, GL_INFO_LOG_LENGTH , &maxLength);

		cerr << prefix
			<< " -- Compiler Log: "
			<< endl
			;

		if (maxLength > 1)
		{
			GLchar * log = new GLchar[maxLength];
			glGetShaderInfoLog(shader, maxLength, 0, log);
			cerr << log
				;
			delete[] log;
		}
	}

	return (compiled != 0);
}
/**
    Returns, whether a shader has been linked with success

    @param	shader		shader id
	@param prefix 		custom error message

*/	
bool Shader::isProgramLinked(GLuint program, char const * prefix)
{
	using namespace std;
	GLint linked(0);

	glGetObjectParameterivARB(program, GL_LINK_STATUS, &linked);

	if (!linked)
	{
		GLint maxLength(0);

		glGetShaderiv(program, GL_INFO_LOG_LENGTH , &maxLength);

		cerr << prefix
			<< " -- Linker Log: "
			<< endl
			;

		if (maxLength > 1)
		{
			GLchar * log = new GLchar[maxLength];
			glGetProgramInfoLog(program, maxLength, 0, log);
			cerr << log
				;
			delete[] log;
		}
	}

	return (linked != 0);
}

bool Shader::shdrSupport = false;

PFNGLCREATESHADERPROC Shader::glCreateShader;
PFNGLCREATEPROGRAMPROC Shader::glCreateProgram;
PFNGLATTACHSHADERPROC Shader::glAttachShader;
PFNGLCOMPILESHADERPROC Shader::glCompileShader;
PFNGLUSEPROGRAMPROC Shader::glUseProgram;
PFNGLDETACHSHADERPROC Shader::glDetachShader;
PFNGLDELETESHADERPROC Shader::glDeleteShader;
PFNGLLINKPROGRAMPROC Shader::glLinkProgram;
PFNGLSHADERSOURCEPROC Shader::glShaderSource;
PFNGLDELETEPROGRAMPROC Shader::glDeleteProgram;
PFNGLGETUNIFORMLOCATIONPROC Shader::glGetUniformLocation;
PFNGLUNIFORM1FPROC Shader::glUniform1f;

PFNGLGETOBJECTPARAMETERIVARBPROC Shader::glGetObjectParameterivARB;
PFNGLGETSHADERIVPROC Shader::glGetShaderiv;
PFNGLGETSHADERINFOLOGPROC Shader::glGetShaderInfoLog;
PFNGLGETPROGRAMINFOLOGPROC Shader::glGetProgramInfoLog;
