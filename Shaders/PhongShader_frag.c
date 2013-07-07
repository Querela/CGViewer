uniform int numLights;
varying vec3 normal, vertex;

void main( )
{
    vec4 color = vec4 ( 0.0, 0.0, 0.0, 1.0 );
    for ( int l = 0 ; l < numLights ; ++ l )
    {
        vec3 lightPosition = gl_LightSource[ l ].position.xyz;
        vec3 vertexToLight = normalize ( lightPosition - vertex );
        color.r += dot ( normal, vertexToLight );
        color.g += dot ( normal, vertexToLight );
        color.b += dot ( normal, vertexToLight );
    }
    
    gl_FragColor = color;
}
