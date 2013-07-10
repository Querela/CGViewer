//#version 150

uniform int numLights;
//uniform sampler2D Texture0;
uniform sampler2D basemap;
uniform sampler2D normalmap;

varying vec3 normal;
varying vec3 vertex;

void main( )
{
    vec4 color = vec4 ( 0.0, 0.0, 0.0, 1.0 );

    vec3 invDir = normalize ( -vertex );

    color += gl_FrontLightModelProduct.sceneColor;
    color += gl_LightModel.ambient * gl_FrontMaterial.ambient;
    //color += gl_LightModel.ambient;

    vec3 N = normal;
//*
    N = vec3 ( texture2D ( normalmap, gl_TexCoord[0].st ) );

    // heuristic to compute tangent and bitangent for normal map
    vec3 Q1 = dFdx ( vertex );
    vec3 Q2 = dFdy ( vertex );
    vec2 st1 = dFdx ( vec2 ( gl_TexCoord[0] ) );
    vec2 st2 = dFdy ( vec2 ( gl_TexCoord[0] ) );
    vec3 tangent = normalize ( Q1 * st2.t - Q2 * st1.t );
    vec3 bitangent = normalize ( -Q1 * st2.s + Q2 * st1.s );

    // smoothing to avoid edges
    tangent = normalize ( cross ( tangent, N ) );
    bitangent = normalize ( cross ( bitangent, N ) );

    // the transpose of texture-to-eye space matrix
    //mat3 TBN = transpose ( mat3 ( tangent, bitangent, normal ) );
    mat3 TBN = mat3 ( tangent, bitangent, normal );

    N = N * TBN; // TBN * N if transpose ?
    //N = normalize ( N );
    N = normalize ( normal + N );
//*/

    int l = 0;
//    for ( int l = 0 ; l < numLights ; ++ l )
//    {
        vec3 lightPosition = gl_LightSource[ l ].position.xyz;
        vec3 vertexToLight = normalize ( lightPosition - vertex );
        float distanceToLight = length ( lightPosition - vertex );

        // too much lightning ?
        //vec4 ambient  = gl_LightSource[l].ambient * gl_FrontMaterial.ambient;

        // different sign if texture ?
        vec3 reflected = reflect ( -vertexToLight, N );

        vec4 diffuse  = gl_LightSource[l].diffuse * gl_FrontMaterial.diffuse * 
                            max ( dot ( N, vertexToLight ), 0.0 );
        vec4 specular = gl_LightSource[l].specular * gl_FrontMaterial.specular *
                            pow ( max ( dot ( reflected, invDir ), 0.0 ), gl_FrontMaterial.shininess );

        float fatt = gl_LightSource[l].constantAttenuation  +
                     gl_LightSource[l].linearAttenuation    * distanceToLight +
                     gl_LightSource[l].quadraticAttenuation * distanceToLight * distanceToLight;
        if ( fatt != 0.0 )
            fatt = 1.0 / fatt;
        else
            fatt = 1.0;

        //color += ambient + diffuse + specular;
        color += fatt * ( diffuse + specular );
//    }

    gl_FragColor = color * texture2D ( basemap, gl_TexCoord[0].st );
}

