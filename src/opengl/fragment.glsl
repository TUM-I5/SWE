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

#define MAX_H 0.5
#define USE_BRIGHTNESS 0

varying vec3 N;
varying vec4 ambient;
varying vec4 worldCoordinates;

/* Shading of water surfaces with fresnel effect */

void main()
{
	// Helper variable for color
	vec4 color = ambient;
	
	// Calculate vector from vertex to the eye, 
	// which resides at (0,0,0)
	vec3 V = -vec3(gl_ModelViewMatrix*worldCoordinates);
	
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
	
	// Scale and limit height
	float h = worldCoordinates.y / 5.0;
	if (h > MAX_H)
		h = MAX_H;
	else if (h < -MAX_H)
		h = -MAX_H;
		
#if (USE_BRIGHTNESS == 1)
	// brightness depending on height
	vec4 blue = vec4(0.2+h, 0.4+h, 0.9+h, 1.0);
#else // USE_BRIGHTNESS
	// color depending on height
	vec4 blue = 2 * vec4(0.5+h, 0.5-abs(h), 0.5-h, abs(h));
	blue = mix(blue, vec4(0.2, 0.4, 0.9, 1.0), 0.3);
#endif // USE_BRIGHTNESS
	color = mix(color, blue, 0.7);
    
    // Set constant transparency 
	gl_FragColor = vec4(color.xyz, 0.7);
}
