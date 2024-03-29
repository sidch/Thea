uniform vec4 material;  // [ka, kl, <ignored>, <ignored>]
uniform vec3 ambient_color;
uniform vec3 light_dir;  // must be specified in camera space, pointing towards object
uniform vec3 light_color;
uniform float two_sided;
uniform float flat_shaded;

varying vec3 position;  // position in camera space
varying vec3 normal;  // normal in camera space

void main()
{
  vec3 N = normalize(flat_shaded > 0.5 ? cross(dFdx(position), dFdy(position)) : normal);
  vec3 L = normalize(light_dir);

  vec3 ambt_color = material[0] * gl_Color.rgb * ambient_color;

  float NdL = -dot(N, L);
  vec3 lamb_color = (NdL >= -two_sided) ? material[1] * abs(NdL) * gl_Color.rgb * light_color : vec3(0.0, 0.0, 0.0);

  gl_FragColor = vec4(ambt_color + lamb_color, 1.0);
}
