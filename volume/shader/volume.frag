#version 150
uniform sampler2D transfer_texture;
uniform sampler3D volume_texture;

uniform float   sampling_distance;
uniform float   sampling_distance_ref;
uniform float   iso_value;
uniform vec3    max_bounds;
uniform ivec3   volume_dimensions;
uniform vec3    obj_to_tex;


uniform vec3    light_ambient_color;
uniform vec3    light_diffuse_color;
uniform vec3    light_specular_color;
uniform float   light_ref_coef;

in		vec3	out_light_ws;

in      vec3      oTexCoord;
in      vec4      oColor;

in per_vertex{
	smooth vec4 ray_exit_os;
	smooth vec4 ray_exit_vs;
	smooth vec4 ray_exit_ts;
} v_in;

in vec4  out_camera_os;
in float out_camera_distance;


out     vec4      FragColor;


struct ray
{
	vec3    origin;
	vec3    direction;
	vec3    direction_rec;
}; // struct ray



/* http://ompf.org/ray/ray_box.html */
//bool_t ray_box_intersects_slabs_geimer_mueller(const aabb_t &box, const ray_t &ray) {
bool
ray_box_intersection(in ray    r,
	in vec3   bbmin,
	in vec3   bbmax,
	out float tmin,
	out float tmax)
{
	float l1 = (bbmin.x - r.origin.x) * r.direction_rec.x;
	float l2 = (bbmax.x - r.origin.x) * r.direction_rec.x;
	tmin = min(l1, l2);
	tmax = max(l1, l2);

	l1 = (bbmin.y - r.origin.y) * r.direction_rec.y;
	l2 = (bbmax.y - r.origin.y) * r.direction_rec.y;
	tmin = max(min(l1, l2), tmin);
	tmax = min(max(l1, l2), tmax);

	l1 = (bbmin.z - r.origin.z) * r.direction_rec.z;
	l2 = (bbmax.z - r.origin.z) * r.direction_rec.z;
	tmin = max(min(l1, l2), tmin);
	tmax = min(max(l1, l2), tmax);

	//return ((lmax > 0.f) & (lmax >= lmin));
	//return ((lmax > 0.f) & (lmax > lmin));
	return ((tmin > 0.0) && (tmax > tmin));
}

bool
inside_volume_bounds(const in vec3 sampling_position)
{
	return (all(greaterThanEqual(sampling_position, vec3(0.0)))
		&& all(lessThanEqual(sampling_position, max_bounds)));
}


float
get_sample_data(vec3 in_sampling_pos) {
	vec3 obj_to_tex = vec3(1.0) / max_bounds;
	return texture(volume_texture, in_sampling_pos * obj_to_tex).r;
}

vec3
get_gradient(vec3 in_sampling_pos, vec3 step_size) {

	
	float nx = get_sample_data(in_sampling_pos - vec3(step_size.x, 0.0, 0.0));
	float px = get_sample_data(in_sampling_pos + vec3(step_size.x, 0.0, 0.0));

	float ny = get_sample_data(in_sampling_pos - vec3(0.0, step_size.z, 0.0));
	float py = get_sample_data(in_sampling_pos + vec3(0.0, step_size.z, 0.0));

	float nz = get_sample_data(in_sampling_pos - vec3(0.0, 0.0, step_size.z));
	float pz = get_sample_data(in_sampling_pos + vec3(0.0, 0.0, step_size.z));

	float grad_x = px - nx;
	float grad_y = py - ny;
	float grad_z = pz - nz;

	return(vec3(grad_x, grad_y, grad_z) * 0.5);

}

vec3
get_gradient_alpha_transfer(vec3 in_sampling_pos, vec3 step_size) {
		
	float nx = get_sample_data(in_sampling_pos - vec3(step_size.x, 0.0, 0.0));
	float px = get_sample_data(in_sampling_pos + vec3(step_size.x, 0.0, 0.0));

	float ny = get_sample_data(in_sampling_pos - vec3(0.0, step_size.z, 0.0));
	float py = get_sample_data(in_sampling_pos + vec3(0.0, step_size.z, 0.0));

	float nz = get_sample_data(in_sampling_pos - vec3(0.0, 0.0, step_size.z));
	float pz = get_sample_data(in_sampling_pos + vec3(0.0, 0.0, step_size.z));

	float dnx = texture(transfer_texture, vec2(nx, 0.5)).a;
	float dpx = texture(transfer_texture, vec2(px, 0.5)).a;

	float dny = texture(transfer_texture, vec2(ny, 0.5)).a;
	float dpy = texture(transfer_texture, vec2(py, 0.5)).a;

	float dnz = texture(transfer_texture, vec2(nz, 0.5)).a;
	float dpz = texture(transfer_texture, vec2(pz, 0.5)).a;

	float grad_x = dpx - dnx;
	float grad_y = dpy - dny;
	float grad_z = dpz - dnz;

	return(vec3(grad_x, grad_y, grad_z) * 0.5);

}


vec4
shade_sample(in vec3 spos, in vec3 grad, in vec4 mat_color, float shadow_mult)
{

	vec3 color = vec3(0.0, 0.0, 0.0);

	if (length(grad) >= 0.0001) {
		
		vec3 lt = out_light_ws;

		vec3 n = normalize(grad);
		vec3 l = normalize(lt - spos);
		vec3 r = reflect(-l, n);

		vec3 v = normalize(spos - out_camera_os.xyz);

		float lambertian = max(dot(l, n), 0.0);

		float specular = 0.0;

		if (lambertian > 0.0) {
			float specAngle = max(dot(r, v), 0.0);
			specular = pow(specAngle, light_ref_coef);
		}
		color.rgb = mat_color.rgb *light_ambient_color +
			mat_color.rgb *lambertian*light_diffuse_color +
			specular*light_specular_color;

	}
	else {
		color.rgb = mat_color.rgb * light_ambient_color;
	}

	color *= shadow_mult;

	return vec4(color, mat_color.a);
}

bool
shadow(vec3 spos) {


	vec3 lt = -out_light_ws;
	vec3 l = normalize(lt - spos);

	vec3 ray_light_increment = l * sampling_distance * 2.0;
	vec3 sampling_pos = spos;
	sampling_pos += ray_light_increment * 3.0;

	float last_data_sample = get_sample_data(sampling_pos);
	float trav_sign = sign(last_data_sample - iso_value);

	bool inside_volume = true;

	while (inside_volume) {

		float cur_data = get_sample_data(sampling_pos);
		float cur_trav_sign = sign(cur_data - iso_value);

		if (cur_trav_sign != trav_sign) {
			return true;
		}

		// increment ray
		sampling_pos += ray_light_increment;
		inside_volume = inside_volume_bounds(sampling_pos);

		last_data_sample = cur_data;
	}

	return false;
}


void main()
{
	ray view_ray;
	view_ray.origin = out_camera_os.xyz;
	view_ray.direction = normalize(v_in.ray_exit_os.xyz - out_camera_os.xyz);
	view_ray.direction_rec = 1.0 / view_ray.direction;
	float t_entry = -1.0;
	float t_exit = -1.0;

	ray_box_intersection(view_ray, vec3(0.0), max_bounds.xyz, t_entry, t_exit);
	if (t_entry < 0.0)
		t_entry = t_exit - length(v_in.ray_exit_vs.xyz);
	
	/// One step trough the volume
	vec3 ray_increment = view_ray.direction * sampling_distance;
	/// Position in Volume
	vec3 sampling_pos = view_ray.origin + view_ray.direction * t_entry + ray_increment; // increment just to be sure we are in the volume

	/// Init color of fragment
	vec4 dst = vec4(0.0, 0.0, 0.0, 0.0);

	/// check if we are inside volume
	bool inside_volume = inside_volume_bounds(sampling_pos);

#if 0
	vec4 max_val = vec4(0.0, 0.0, 0.0, 0.0);

	/// the traversal loop,
	/// termination when the sampling position is outside volume boundarys
	/// another termination condition for early ray termination is added
	while (inside_volume)
	{
		/// get sample
		float s = get_sample_data(sampling_pos);

		/// apply the transfer functions to retrieve color and opacity
		vec4 color = texture(transfer_texture, vec2(s, s));
		
		max_val = max(vec4(s), max_val);
		
		/// increment the ray sampling position
		sampling_pos += ray_increment;

		/// update the loop termination condition
		inside_volume = inside_volume_bounds(sampling_pos);		
	}

	dst = max_val;

#else

	float last_data_sample = get_sample_data(sampling_pos);
	
	float trav_sign = sign(last_data_sample - iso_value);
	float cur_trav_sign = trav_sign;
	
	bool found_intersection = false;

	while (inside_volume && !found_intersection) {

		
		float cur_data = get_sample_data(sampling_pos);
		trav_sign = cur_trav_sign;
		float cur_trav_sign = sign(cur_data - iso_value);

		if (cur_trav_sign != trav_sign) {
			// ok a sign change has occured
			// we crossed the iso surface
			while (length(ray_increment) > 0.0001) {
				// first try:
				//  - interpolate the intrersection position
				//float u = (iso_value - cur_data) / (last_data_sample - cur_data);
				//vec3 iso_intersect = mix((sampling_pos - ray_increment), sampling_pos, u);

				float cur_data = get_sample_data(sampling_pos);
				cur_trav_sign = sign(cur_data - iso_value);
				
				if (trav_sign != cur_trav_sign) {
					ray_increment *= -1.0;
				}
				ray_increment *= 0.5;

				trav_sign = cur_trav_sign;

				sampling_pos += ray_increment;
			}

			//float u = (iso_value - cur_data) / (last_data_sample - cur_data);
			//vec3 iso_intersect = mix((sampling_pos - ray_increment), sampling_pos, u);
			vec3 iso_intersect = sampling_pos;// mix((sampling_pos - ray_increment), sampling_pos, u);
			//dst.rgb = iso_intersect;
			//dst.a = 1.0;

			float shadow_mult = 1.0;
			if (shadow(sampling_pos)) {
				shadow_mult = 0.2;
			}
			vec3 gradient = get_gradient(iso_intersect, vec3(0.01, 0.01, 0.01));
			vec4 src = shade_sample(iso_intersect, gradient, vec4(1.0, 1.0, 1.0, 1.0), shadow_mult);

			dst = src;

			found_intersection = true;

		}
		// increment ray
		sampling_pos += ray_increment;
		inside_volume = inside_volume_bounds(sampling_pos);

		last_data_sample = cur_data;
	}
#endif
	FragColor = dst;
	//FragColor = oColor * texture2D(transfer_texture, oTexCoord.xy);	
	//FragColor = vec4(ray_entry_position, 0.8);
	//FragColor = vec4(abs(camera_location), 0.8);
	//FragColor = vec4(normalize(abs(ray_increment)), 1.0);
	//FragColor.rgb = out_camera_os.rgb;
	//FragColor.a = 1.0;
	return;
};
