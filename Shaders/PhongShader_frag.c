//#version 150

uniform int numLights;
uniform sampler2D basemap;
uniform sampler2D normalmap;

varying vec3 normal;
varying vec3 vertex;

// constants for toon effects
const float outline_threshold = 0.2;

void main( )
{
    vec4 color = vec4 ( 0.0, 0.0, 0.0, 1.0 );
    color += gl_FrontLightModelProduct.sceneColor;
    color += gl_LightModel.ambient * gl_FrontMaterial.ambient;

    vec3 invDir = normalize ( -vertex );

    vec3 N = normal;
//*
    N = vec3 ( texture2D ( normalmap, gl_TexCoord[0].st ) );

    // heuristic to compute tangent and bitangent for normal map
    vec3 Q1 = dFdx ( vertex );
    vec3 Q2 = dFdy ( vertex );
    vec2 st1 = dFdx ( gl_TexCoord[0].st );
    vec2 st2 = dFdy ( gl_TexCoord[0].st );
    vec3 tangent = normalize ( Q1 * st2.t - Q2 * st1.t );
    vec3 bitangent = normalize ( -Q1 * st2.s + Q2 * st1.s );

    // smoothing to avoid edges
    tangent = normalize ( cross ( tangent, N ) );
    bitangent = normalize ( cross ( bitangent, N ) );

    // the transpose of texture-to-eye space matrix
    //mat3 TBN = transpose ( mat3 ( tangent, bitangent, normal ) );
    mat3 TBN = mat3 ( tangent, bitangent, normal );

    N = N * TBN; // TBN * N if transpose ?
    N = normalize ( normal + N );
//*/

    vec4 texture_color = texture2D ( basemap, gl_TexCoord[0].st );

    int l = 0;
//    for ( int l = 0 ; l < numLights ; ++ l )
//    {
        vec3 lightPosition = gl_LightSource[ l ].position.xyz;
        vec3 vertexToLight = normalize ( lightPosition - vertex );
        float distanceToLight = length ( lightPosition - vertex );

        vec3 reflected = reflect ( -vertexToLight, N );

        // use intensity of diffuse light to create light steps
        float I_diffuse =  max ( dot ( N, vertexToLight ), 0.0 );
        if      ( I_diffuse > 0.8 ) color = vec4 ( 1.0 );
        else if ( I_diffuse > 0.7 ) color = vec4 ( vec3 ( 0.6 ), 1.0 );
        else if ( I_diffuse > 0.4 ) color = vec4 ( vec3 ( 0.2 ), 1.0 );
        else if ( I_diffuse > 0.1 ) color = vec4 ( vec3 ( 0.05 ), 1.0 );

//        vec4 specular = gl_LightSource[l].specular * gl_FrontMaterial.specular *
//                            pow ( max ( dot ( reflected, invDir ), 0.0 ), gl_FrontMaterial.shininess );
//
//        float fatt = gl_LightSource[l].constantAttenuation  +
//                     gl_LightSource[l].linearAttenuation    * distanceToLight +
//                     gl_LightSource[l].quadraticAttenuation * distanceToLight * distanceToLight;
//        if ( fatt != 0.0 )
//            fatt = 1.0 / fatt;
//        else
//            fatt = 1.0;

//    }

    // use the texture color to get some other shades than black and white
    gl_FragColor = color * smoothstep ( 0.3, 0.6, texture_color );

    // outline
    if ( dot ( invDir, N ) < outline_threshold )
        gl_FragColor = vec4 ( 0.0, 0.0, 0.0, 1.0 );
}

