varying vec3 position;  // position in camera space
varying vec3 normal;  // normal in camera space

void main()
{
  gl_Position = ftransform();

  position = vec3(gl_ModelViewMatrix * gl_Vertex);  // assume rigid transform, so we can drop w
  normal = gl_NormalMatrix * gl_Normal;

  gl_FrontColor = gl_Color;
  gl_BackColor = gl_Color;
}
