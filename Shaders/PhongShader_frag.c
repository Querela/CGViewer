uniform int numLights;
uniform sampler2D Texture0;

varying vec3 normal;
varying vec3 vertex;

void main( )
{
    //vec4 color = vec4 ( 0.0, 0.0, 0.0, 1.0 );
/*
    for ( int l = 0 ; l < numLights ; ++ l )
    {
        vec3 lightPosition = gl_LightSource[ l ].position.xyz;
        vec3 vertexToLight = normalize ( lightPosition - vertex );

        color.r += dot ( normal, vertexToLight );
        color.g += dot ( normal, vertexToLight );
        color.b += dot ( normal, vertexToLight );
    }
*/
    gl_FragColor = texture2D ( Texture0, vec2 ( gl_TexCoord[0] ) );
}
