uniform int numLights;
uniform sampler2D Texture0;
uniform sampler2D basemap;

varying vec3 normal;
varying vec3 vertex;

void main( )
{
    vec4 color = vec4 ( 0.0, 0.0, 0.0, 1.0 );

    vec3 invDir = normalize ( -vertex );

    //color += gl_FrontLightModelProduct.sceneColor;
    //color += gl_LightModel.ambient * gl_FrontMaterial.ambient;
    color += gl_LightModel.ambient;

    int l = 0;
//    for ( int l = 0 ; l < numLights ; ++ l )
//    {
        vec3 lightPosition = gl_LightSource[ l ].position.xyz;
        vec3 vertexToLight = normalize ( lightPosition - vertex );
        float distanceToLight = length ( lightPosition - vertex );

        vec3 reflected = normalize ( reflect ( -vertexToLight, normal ) );

        vec4 ambient  = gl_LightSource[l].ambient * gl_FrontMaterial.ambient;
        vec4 diffuse  = gl_LightSource[l].diffuse * gl_FrontMaterial.diffuse * 
                            max ( dot (normal, vertexToLight ), 0.0 );
        vec4 specular = gl_LightSource[l].specular * gl_FrontMaterial.specular *
                            pow ( max ( dot ( reflected, invDir ), 0.0 ), gl_FrontMaterial.shininess );

        float fatt = 1.0;
        if ( gl_LightSource[l].constantAttenuation != 0.0 &&
             gl_LightSource[l].linearAttenuation != 0.0 &&
             gl_LightSource[l].quadraticAttenuation != 0.0 )
        {
            fatt = 1.0 / ( gl_LightSource[l].constantAttenuation +
                           gl_LightSource[l].linearAttenuation * distanceToLight +
                           gl_LightSource[l].quadraticAttenuation * distanceToLight * distanceToLight );
        }

        //color += ambient + diffuse + specular;
        color += fatt * ( diffuse + specular );
//    }

    gl_FragColor = color * texture2D ( Texture0, vec2 ( gl_TexCoord[0] ) );
}
