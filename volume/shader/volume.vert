#version 150
#extension GL_ARB_shading_language_420pack : require
#extension GL_ARB_explicit_attrib_location : require

uniform mat4 matWVP;
in      vec4 Position;
in      vec4 Color;
in      vec3 TexCoord;
out     vec3 oTexCoord;
out     vec4 oColor;

uniform mat4    Model;
uniform mat4    Modelview;
uniform mat4    Modelview_inv;

uniform vec3 obj_to_tex;

uniform vec3    light_position;

out per_vertex{
	smooth vec4 ray_exit_os;
	smooth vec4 ray_exit_vs;
	smooth vec4 ray_exit_ts;
} vertex_out;

out vec4  out_camera_os;
out float out_camera_distance;
out vec3 out_light_ws;

void main()
{
	gl_Position = (matWVP * Position);
	oTexCoord = TexCoord;
	oColor.rgb = pow(Color.rgb, vec3(2.2));   // convert from sRGB to linear
	oColor.a = Color.a;

	vertex_out.ray_exit_os = vec4(Position.xyz, 1.0);
	vertex_out.ray_exit_vs = Modelview * vec4(Position.xyz, 1.0);
	vertex_out.ray_exit_ts = vec4(Position.xyz, 1.0) * vec4(obj_to_tex, 1.0);
	
	out_camera_os		= Modelview_inv * vec4(0.0, 0.0, 0.0, 1.0);
	out_camera_distance = length(vertex_out.ray_exit_vs.xyz);

	out_light_ws = (Model *  vec4(light_position, 1.0)).xyz;
};
