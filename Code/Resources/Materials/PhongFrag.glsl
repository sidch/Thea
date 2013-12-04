uniform vec3 ambient_color;
uniform vec3 light_dir;  // must be specified in camera space, pointing towards object
uniform vec3 light_color;
uniform vec4 material;  // [ka, kl, <ignored>, <ignored>]
uniform float two_sided;

varying vec3 position;  // position in camera space
varying vec3 normal;  // normal in camera space

void main()
{
  vec3 N = normalize(normal);
  vec3 L = normalize(light_dir);

  vec3 ambt_color = material[0] * gl_Color.rgb * ambient_color;
  // ambt_color = vec3(0.0, 0.0, 0.0);

  float NdL = -dot(N, L);
  vec3 lamb_color = (NdL >= -two_sided) ? material[1] * abs(NdL) * gl_Color.rgb * light_color : vec3(0.0, 0.0, 0.0);
  // lamb_color = vec3(0.0, 0.0, 0.0);

  gl_FragColor = vec4(ambt_color + lamb_color, 1.0);
}
