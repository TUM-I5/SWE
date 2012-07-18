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

varying vec3 V;
varying vec3 N;
varying vec4 ambient;

/* Shading of water surfaces with fresnel effect */

void main()
{
	// Helper variable for color
	vec4 color = ambient;
	
	// Get normal and view vector
	vec3 view = normalize(V);
	vec3 normal = normalize(N); ;
	if (!gl_FrontFacing) {
		// Handle backface rendering
		normal = -normal;
	}
	
	// Compute fresnel factor
	float fresnelFactor = 0.1+0.1*pow(1+dot(normal,view), 3.0);
		
	// Combine color with fresnel factor
	color += gl_LightSource[0].specular * fresnelFactor;
	color = mix(color, vec4(0.2, 0.4, 0.9, 1.0), 0.2);
    
    // Set constant transparency 
	gl_FragColor = vec4(color.xyz, 0.7);
}
