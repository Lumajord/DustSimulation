/*
 * Copyright 1993-2009 NVIDIA Corporation.  All rights reserved.
 *
 * NVIDIA Corporation and its licensors retain all intellectual property and 
 * proprietary rights in and to this software and related documentation and 
 * any modifications thereto.  Any use, reproduction, disclosure, or distribution 
 * of this software and related documentation without an express license 
 * agreement from NVIDIA Corporation is strictly prohibited.
 * 
 */
 
 #define STRINGIFY(A) #A

// vertex shader
const char *vertexShader = STRINGIFY(
uniform float pointRadius;  // point size in world space
uniform float pointScale;   // scale to calculate size in pixels

void main()
{
    // calculate window-space point size
    vec3 posEye = vec3(gl_ModelViewMatrix * vec4(gl_Vertex.xyz, 1.0));
    float dist = length(posEye);
    gl_PointSize = pointRadius * (pointScale / dist);

    gl_TexCoord[0] = gl_MultiTexCoord0;
    gl_Position = gl_ModelViewProjectionMatrix * vec4(gl_Vertex.xyz, 1.0);

    gl_FrontColor = gl_Color;
}
);

// pixel shader for rendering points as shaded spheres
const char *sphereBoundariesPixelShader = STRINGIFY(

uniform vec4 clear_color;
uniform vec4 border_color;

uniform float min_border_depth;
uniform float distance_factor;

void main()
{
    const vec3 lightDir = vec3(0.577, 0.577, 0.577);

    // calculate normal from texture coordinates
    vec3 N;
    N.xy = gl_TexCoord[0].xy*vec2(4.0, -4.0) + vec2(-2, 2);
    float mag = dot(N.xy, N.xy);

    if(mag > 4)
		discard;   // kill pixels outside circle
	else if(mag > 1.0)
	{
		mag = sqrt(mag) - 1.0f;
		
		gl_FragDepth = gl_DepthRange.far / gl_DepthRange.diff  - gl_DepthRange.far * gl_DepthRange.near / gl_DepthRange.diff / ( gl_FragCoord.z/gl_FragCoord.w + min_border_depth + mag * distance_factor);
		gl_FragColor = border_color;
		//float blend_factor = (1.0f - mag * mag);
		//blend_factor * border_color + (1.0f - blend_factor) * clear_color;
	}
	else
	{
		N.z = sqrt(1.0-mag);
		
		// calculate lighting
		float diffuse = max(0.2, dot(lightDir, N));

		gl_FragColor = gl_Color * diffuse;
		gl_FragDepth = gl_FragCoord.z;
	}
}
);

const char *boundariesPixelShader = STRINGIFY(

uniform vec4 clear_color;
uniform vec4 border_color;

uniform float min_border_depth;
uniform float distance_factor;

void main()
{
    const vec3 lightDir = vec3(0.577, 0.577, 0.577);

    // calculate normal from texture coordinates
    vec3 N;
    N.xy = gl_TexCoord[0].xy*vec2(4.0, -4.0) + vec2(-2.0, 2.0);
    float mag = dot(N.xy, N.xy);

    if (mag > 4)
		discard;   // kill pixels outside circle
	else if(mag > 1.0)
	{
		mag = sqrt(mag) - 1.0f;

		// blend border color with clear color
		//float blend_factor = (1.0f - mag * mag);
		//gl_FragColor = blend_factor * border_color + (1.0f - blend_factor) * clear_color;
		gl_FragColor = border_color;

		gl_FragDepth = gl_DepthRange.far / gl_DepthRange.diff  - gl_DepthRange.far * gl_DepthRange.near / gl_DepthRange.diff / ( gl_FragCoord.z/gl_FragCoord.w + min_border_depth + mag * distance_factor);
	}
	else
	{
		gl_FragColor = clear_color;
		gl_FragDepth = gl_FragCoord.z;
	}
}
);

// pixel shader for rendering points as shaded spheres
const char *spherePixelShader = STRINGIFY(

uniform sampler2D noise_texture;

void main()
{
    const vec3 lightDir = vec3(0.577, 0.577, 0.577);

    // calculate normal from texture coordinates
    vec3 N;
    N.xy = gl_TexCoord[0].xy*vec2(2.0, -2.0) + vec2(-1.0, 1.0);
    float mag = dot(N.xy, N.xy);

    if (mag > 1.0)
		discard;   // kill pixels outside circle
	else
	{
		N.z = sqrt(1.0-mag);

		// calculate lighting
		float diffuse = max(0.2, dot(lightDir, N));

		// apply fogging
		//diffuse *= gl_FragCoord.z/gl_FragCoord.z;
		//float fog = zFar / ( zFar - zNear ) * gl_FragCoord.z;
		//float noise = texture2D(noise_texture, gl_TexCoord[0].xy);
		//gl_FragColor = texture2D(noise_texture, gl_TexCoord[0].xy);
		//gl_FragColor = vec4(gl_TexCoord[0].xy, 1, 1);

		gl_FragColor = gl_Color * diffuse;
		//gl_FragDepth = gl_FragCoord.z;
	}
}
);
