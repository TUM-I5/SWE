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

varying vec3 N;
varying vec3 V;
varying vec4 ambient;

/* Shading of water surfaces with fresnel effect */

void main()
{	
	// Compute properties of the light source
	ambient = gl_FrontMaterial.ambient * gl_LightSource[0].ambient;
	ambient += gl_LightModel.ambient * gl_FrontMaterial.ambient;
	
	// Transform our normal into eye space
	N = normalize(gl_NormalMatrix * gl_Normal);
	
	// Calculate vector from vertex to the eye, 
	// which resides at (0,0,0)
	V = -vec3(gl_ModelViewMatrix*gl_Vertex);
	
	// Compute vertex position via internal transform function
	gl_Position = ftransform();
} 