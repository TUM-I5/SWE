// =====================================================================
// This file is part of SWE_CUDA (see file SWE_Block.cu for details).
//
// Copyright (C) 2012 Sebastian Rettenberger
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

#ifndef TEXT_H
#define TEXT_H

#include <cmath>
#include <string>
#include <vector>
#include <SDL/SDL_ttf.h>
#include <SDL/SDL_opengl.h>

class Text
{
private:
	/** OpenGL textures containing the text */
	std::vector<GLuint> textures;
	/** Width of each text */
	std::vector<int> width;
	/** Height of each text */
	std::vector<int> height;

	/** Index of the next text, that should be rendered */
	unsigned int nextText;

public:
	Text()
	{
		// Initialize TTF if not done yet
		if (instances == 0) {
			if (TTF_Init() < 0) {
		        fprintf(stderr, "\nCould not initialize SDL_ttf\n" );
		        exit(-1);
			}

			openFont();
		}

		instances++;
	}

	~Text()
	{
		// We assume that std::vector is just an array ...
		glDeleteTextures(textures.size(), &textures[0]);

		instances--;

		if (instances == 0) {
			TTF_CloseFont(font);
			TTF_Quit();
		}
	}

	void addText(const char* text)
	{
		if (!font)
			return;

		// Draw string
		SDL_Color black = {0, 0, 0, 255};
		SDL_Surface* surf = TTF_RenderText_Blended(font, text, black);

		if (surf) {
			GLuint texture;

			/* Tell GL about our new texture */
			glGenTextures(1, &texture);
			glBindTexture(GL_TEXTURE_2D, texture);
			glTexImage2D(GL_TEXTURE_2D, 0, 4, surf->w, surf->h, 0, GL_BGRA,
					GL_UNSIGNED_BYTE, surf->pixels);

			width.push_back(surf->w);
			height.push_back(surf->h);

			textures.push_back(texture);

			/* Clean up */
			SDL_FreeSurface(surf);
		}
	}

	void startTextMode()
	{
		// Enable 2D mode
		int vPort[4];

		glGetIntegerv(GL_VIEWPORT, vPort);

		glMatrixMode(GL_PROJECTION);
		glPushMatrix();
		glLoadIdentity();

		glOrtho(0, vPort[2], 0, vPort[3], -1, 1);
		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		glLoadIdentity();
		glDisable(GL_DEPTH_TEST);
	    glEnable(GL_BLEND);
		glEnable(GL_TEXTURE_2D);
	    glBlendFunc(GL_ONE_MINUS_DST_COLOR, GL_ONE_MINUS_SRC_ALPHA);

	    nextText = 0;
	}

	/**
	 * @return True there are more textures
	 */
	bool showNextText(SDL_Rect &location)
	{
		if (!font)
			return false;

		if (nextText >= textures.size())
			// No more textures to show
			return false;

		/* GL_NEAREST looks horrible, if scaled... */
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

		/* prepare to render our texture */
		glBindTexture(GL_TEXTURE_2D, textures[nextText]);

		/* Draw a quad at location */
		glBegin(GL_QUADS);
			/* Recall that the origin is in the lower-left corner
			   That is why the TexCoords specify different corners
			   than the Vertex coords seem to. */
			glTexCoord2f(0.0f, 1.0f);
				glVertex2f(location.x, location.y - height[nextText]);
			glTexCoord2f(1.0f, 1.0f);
				glVertex2f(location.x + width[nextText], location.y - height[nextText]);
			glTexCoord2f(1.0f, 0.0f);
				glVertex2f(location.x + width[nextText], location.y);
			glTexCoord2f(0.0f, 0.0f);
				glVertex2f(location.x, location.y);
		glEnd();

		/* return the deltas in the unused w,h part of the rect */
		location.w = width[nextText];
		location.h = height[nextText];

		nextText++;

		return nextText < textures.size();
	}

	void endTextMode()
	{
		// Disable 2D mode
		glDisable(GL_BLEND);
		glDisable(GL_TEXTURE_2D);
		glEnable(GL_DEPTH_TEST);
		glMatrixMode(GL_PROJECTION);
		glPopMatrix();
		glMatrixMode(GL_MODELVIEW);
		glPopMatrix();
	}

private:
	static void openFont()
	{
		std::string fullFileName = "";

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
			path = path.substr(0, t);
		}

		fullFileName = path + std::string("/");
#endif
		fullFileName.append("FreeSans.ttf");

		font = TTF_OpenFont(fullFileName.c_str(), 16);

		if (!font) {
			fprintf(stderr, "TTF_OpenFont: %s\n", TTF_GetError());
			fprintf(stderr, "All text will be ignored\n");
		}
	}

	/** Number of text classes */
	static unsigned int instances;

	/** The only font we use */
	static TTF_Font* font;
};

#endif // TEXT_H
