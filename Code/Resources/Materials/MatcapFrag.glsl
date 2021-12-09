uniform sampler2D matcap_tex;
uniform float two_sided;
uniform float flat_shaded;

varying vec3 position;  // position in camera space
varying vec3 normal;  // normal in camera space

float lightness(vec3 color)
{
  return 0.5 * (max(max(color.r, color.g), color.b) + min(min(color.r, color.g), color.b));
}

void main()
{
  vec3 N = normalize(flat_shaded > 0.5 ? cross(dFdx(position), dFdy(position)) : normal);
  vec2 matcap_uv = (two_sided < 0.5 && N.z < 0.0 ? normalize(N.xy) : N.xy);
  vec4 matcap_color = texture2D(matcap_tex, 0.495 * matcap_uv + 0.5);

  float lite = lightness(matcap_color.rgb);
  float highlight_blend = lite * lite * lite * lite;
  gl_FragColor = vec4(mix(gl_Color.rgb, vec3(1.0, 1.0, 1.0), highlight_blend) * matcap_color.rgb, gl_Color.a);
}
