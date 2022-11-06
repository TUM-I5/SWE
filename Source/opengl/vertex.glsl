// =====================================================================
// This file is part of SWE_CUDA (see file SWE_Block.cu for details).
// 
// Copyright (C) 2010,2011 Tobias Schnabel
// Copyright (C) 2012      Sebastian Rettenberger
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

uniform float scale = 10;
varying vec3 N;
varying vec4 ambient;
varying vec4 worldCoordinates;

/* Shading of water surfaces with fresnel effect */

void main()
{	
	// Compute properties of the light source
	ambient = gl_FrontMaterial.ambient * gl_LightSource[0].ambient;
	ambient += gl_LightModel.ambient * gl_FrontMaterial.ambient;
	
	// Transform our normal into eye space
	N = normalize(gl_NormalMatrix * gl_Normal);
	
	// Save world coordinates	
	worldCoordinates = gl_Vertex * scale;
	
	// Compute vertex position via internal transform function
	gl_Position = ftransform();
} 