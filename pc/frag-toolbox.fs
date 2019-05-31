//Library of helpful fragment shader functions taken from various places.
//Book of Shaders
// https://thebookofshaders.com/
//IQ website
// https://www.iquilezles.org/www/index.htm
//Shader toy compatible but easily transferred to other shaders.

//USEFUL MATH FUNCITONS: inigo
float almostIdentity( float x, float m, float n )
{
    if( x>m ) return x;

    const float a = 2.0*n - m
    const float b = 2.0*m - 3.0*n;
    const float t = x/m;

    return (a*t + b)*t*t + n;
}

float impulse( float k, float x )
{
    const float h = k*x;
    return h*exp(1.0-h);
}

float cubicPulse( float c, float w, float x )
{
    x = fabs(x - c);
    if( x>w ) return 0.0;
    x /= w;
    return 1.0 - x*x*(3.0-2.0*x);
}


float expStep( float x, float k, float n )
{
    return exp( -k*pow(x,n) );
}


float gain(float x, float k) 
{
    float a = 0.5*pow(2.0*((x<0.5)?x:1.0-x), k);
    return (x<0.5)?a:1.0-a;
}


float parabola( float x, float k )
{
    return pow( 4.0*x*(1.0-x), k );
}


float pcurve( float x, float a, float b )
{
    float k = pow(a+b,a+b) / (pow(a,a)*pow(b,b));
    return k * pow( x, a ) * pow( 1.0-x, b );
}

//bounce
float sinc( float x, float k )
{
    const float a = PI * ((float(k)*x-1.0);
    return sin(a)/a;
}

// Color
//Red Green Blue -> Hue Saturation Brightness
vec3 rgb2hsb( in vec3 c ){
    vec4 K = vec4(0.0, -1.0 / 3.0, 2.0 / 3.0, -1.0);
    vec4 p = mix(vec4(c.bg, K.wz),
                 vec4(c.gb, K.xy),
                 step(c.b, c.g));
    vec4 q = mix(vec4(p.xyw, c.r),
                 vec4(c.r, p.yzx),
                 step(p.x, c.r));
    float d = q.x - min(q.w, q.y);
    float e = 1.0e-10;
    return vec3(abs(q.z + (q.w - q.y) / (6.0 * d + e)),
                d / (q.x + e),
                q.x);
}

//  Function from Iñigo Quiles
//  https://www.shadertoy.com/view/MsS3Wc
vec3 hsb2rgb( in vec3 c ){
    vec3 rgb = clamp(abs(mod(c.x*6.0+vec3(0.0,4.0,2.0),
                             6.0)-3.0)-1.0,
                     0.0,
                     1.0 );
    rgb = rgb*rgb*(3.0-2.0*rgb);
    return c.z * mix(vec3(1.0), rgb, c.y);
}

//noise 
// Based on Morgan McGuire @morgan3d
// https://www.shadertoy.com/view/4dS3Wd
float noise (in vec2 st) {
    vec2 i = floor(st);
    vec2 f = fract(st);

    // Four corners in 2D of a tile
    float a = random(i);
    float b = random(i + vec2(1.0, 0.0));
    float c = random(i + vec2(0.0, 1.0));
    float d = random(i + vec2(1.0, 1.0));

    vec2 u = f * f * (3.0 - 2.0 * f);

    return mix(a, b, u.x) +
            (c - a)* u.y * (1.0 - u.x) +
            (d - b) * u.x * u.y;
}


//Structure
//ro: camera position, rd: camera direction
vec3 render(vec3 ro, vec rd){
	//computer color based on raytracing/whatever here
	return col;
}

//Anti-Aliasing
// AA
for( int m=ZERO; m<AA; m++ )
    for( int n=ZERO; n<AA; n++ )
    {
		//compute total color
		//tot += col;
	}
    tot /= float(AA*AA);
	}
	
//CAMERAS
// Look At Camera
//ro: camera position, ta: look at point, cr: fov angle
mat3 setCamera( in vec3 ro, in vec3 ta, float cr )
{
	vec3 cw = normalize(ta-ro);
	vec3 cp = vec3(sin(cr), cos(cr),0.0);
	vec3 cu = normalize( cross(cw,cp) );
	vec3 cv =          ( cross(cu,cw) );
    return mat3( cu, cv, cw );				//rotational camera basis
}

//Alternative Camera
    vec2 m = vec2(0.5);
	float an = sin(-0.25 + 0.31416*iTime) - 6.2831*(m.x-0.5);

	vec3 ro = vec3(3.5*sin(an),1.8,3.5*cos(an));
    vec3 ta = vec3(0.0,1.5,0.0);

    // camera matrix
    vec3 ww = normalize( ta - ro );
    vec3 uu = normalize( cross(ww,vec3(0.0,1.0,0.0) ) );
    vec3 vv = normalize( cross(uu,ww));

	// create view ray
	vec3 rd = normalize( p.x*uu + p.y*vv + 2.0*ww );

//Screen Coordinates
vec2 q = fragCoord.xy / iResolution.xy;
vec2 p = -1.0 + 2.0 * q;
p.x *= iResolution.x/iResolution.y;

//VISUALIZATION
// Plot a line on Y using a value between 0.0-1.0
float plot(vec2 st, float pct){
  return  smoothstep( pct-0.1, pct, st.y) -
          smoothstep( pct, pct+0.1, st.y);
}


//COMPUTE CYLINDER AND FOG

float computeFog( vec3  ro, vec3  rd,   // ray origin, ray direction
                  vec3  sc, float sr,   // sphere center, sphere radius
                  float dbuffer )
{
    // normalize the problem to the canonical sphere
    float ndbuffer = dbuffer / sr;
    vec3  rc = (ro - sc)/sr;
	
    // find intersection with sphere
    float b = dot(rd,rc);
    float c = dot(rc,rc) - 1.0f;
    float h = b*b - c;

    // not intersecting
    if( h<0.0f ) return 0.0f;
	
    h = sqrtf( h );
    float t1 = -b - h;
    float t2 = -b + h;

    // not visible (behind camera or behind ndbuffer)
    if( t2<0.0f || t1>ndbuffer ) return 0.0f;

    // clip integration segment from camera to ndbuffer
    t1 = max( t1, 0.0f );
    t2 = min( t2, ndbuffer );

    // analytical integration of an inverse squared density
    float i1 = -(c*t1 + b*t1*t1 + t1*t1*t1/3.0f);
    float i2 = -(c*t2 + b*t2*t2 + t2*t2*t2/3.0f);
    return (i2-i1)*(3.0f/4.0f);
}

float computeFogCylinder( vec3  cd, vec3  rd,   // camera direction, ray direction
                  vec3  r2c, float cr,   // point to camera, cylinder radius
                  float depth )
{
// Cylinder intersection
//ro = camera position
//rd = camera vector
//co = cylinder position
//cr = cylinder radius
//depth = object depth at pixel
//cd = cylinder orientation
//r2c = nearest projected point ALONG cylinder axis
//          (dot((ro-co),cd)*cd+co)-co

float3 rayproj = rd - dot(rd, cd) * cd;
//float3 rayproj = float3(rd.x,rd.y,0.0);

float a = dot(rayproj, rayproj);
float b = dot(r2c, rayproj);
float c = dot(r2c, r2c) - cr * cr;

float density = 0.0;


// Find intersection
float h = b*b - a * c;


h = max(h,0);
h = sqrt(h);

//enter
float t1 = -(h+b) /a;
//exit
float t2 = (h-b) /a;

t1 = max(t1, 0);
t2 = min(t2, depth);

float bValid = (t2 >= 0.0 && t1 <= depth);


// Inverse Squared Density
float i1 = -(c*t1 + b*t1*t1 + a*t1*t1*t1/3.0);
float i2 = -(c*t2 + b*t2*t2 + a*t2*t2*t2/3.0);

density = (i2-i1) / (cr * cr);
density = min(density, 200.f);

density *= bValid;

float scale =  (3.0 / 4.0) / cr;


return density * scale;

}

float computeFogCone( vec3  cd, vec3  rd,   // camera direction, ray direction
                  vec3  r2c, float cr,   // point to camera, cylinder radius
                  float depth )
{
// Cylinder intersection
//ro = camera position
//rd = camera vector
//co = cylinder position
//cr = cylinder radius
//depth = object depth at pixel
//cd = cylinder orientation
//r2c = nearest projected point ALONG cylinder axis
//          (dot((ro-co),cd)*cd+co)-co

float3 rayproj = rd - dot(rd, cd) * cd;
//float3 rayproj = float3(rd.x,rd.y,0.0);

float a = dot(rayproj, rayproj);
float b = dot(r2c, rayproj);
float c = dot(r2c, r2c) - cr * cr;

float density = 0.0;


// Find intersection
float h = b*b - a * c;


h = max(h,0);
h = sqrt(h);

//enter
float t1 = -(h+b) /a;
//exit
float t2 = (h-b) /a;

t1 = max(t1, 0);
t2 = min(t2, depth);

float bValid = (t2 >= 0.0 && t1 <= depth);


// Inverse Squared Density
float i1 = -(c*t1 + b*t1*t1 + a*t1*t1*t1/3.0);
float i2 = -(c*t2 + b*t2*t2 + a*t2*t2*t2/3.0);

density = (i2-i1) / (cr * cr);
density = min(density, 200.f);

density *= bValid;

float scale =  (3.0 / 4.0) / cr;


return density * scale;

}

//https://www.shadertoy.com/view/MtcXWr
//https://www.shadertoy.com/view/XsXSzj
//http://lousodrome.net/blog/light/2017/01/03/intersection-of-a-ray-and-a-cone/

// cone inscribed in a unit cube centered at 0
bool cone(vec3 org, vec3 dir, out float near, out float far)
{
	// scale and offset into a unit cube
	org.x += 0.5;
	float s = 0.5;
	org.x *= s;
	dir.x *= s;
	
	// quadratic x^2 = y^2 + z^2
	float a = dir.y * dir.y + dir.z * dir.z - dir.x * dir.x;
	float b = org.y * dir.y + org.z * dir.z - org.x * dir.x;
	float c = org.y * org.y + org.z * org.z - org.x * org.x;
    
	float cap = (s - org.x) / dir.x;
	
	// linear
	if( a == 0.0 )
	{
		near = -0.5 * c/b;
		float x = org.x + near * dir.x;
		if( x < 0.0 || x > s )
			return false; 

		far = cap;
		float temp = min(far, near); 
		far = max(far, near);
		near = temp;
		return far > 0.0;
	}

	float delta = b * b - a * c;
	if( delta < 0.0 )
		return false;

	// 2 roots
	float deltasqrt = sqrt(delta);
	float arcp = 1.0 / a;
	near = (-b - deltasqrt) * arcp;
	far = (-b + deltasqrt) * arcp;
	
	// order roots
	float temp = min(far, near);
	far = max(far, near);
	near = temp;

	float xnear = org.x + near * dir.x;
	float xfar = org.x + far * dir.x;

	if( xnear < 0.0 )
	{
		if( xfar < 0.0 || xfar > s )
			return false;
		
		near = far;
		far = cap;
	}
	else if( xnear > s )
	{
		if( xfar < 0.0 || xfar > s )
			return false;
		
		near = cap;
	}
	else if( xfar < 0.0 )
	{
		// The apex is problematic,
		// additional checks needed to
		// get rid of the blinking tip here.
		far = near;
		near = cap;
	}
	else if( xfar > s )
	{
		far = cap;
	}
	
	return far > 0.0;
}

//-------------------------------- ---------------------------------------------
//----------------------------CONE----------------------------------------------
//-------------------------------- ---------------------------------------------

// cone inscribed in a unit cube centered at 0
bool cone(vec3 org, vec3 dir, out float near, out float far, out float density)
{
	// scale and offset into a unit cube
	//org.x += 0.5;
	float s = 0.5;
	//org.x *= s;
	//dir.x *= s;
	
    
    vec3 center = vec3(0.0,0.0,0.0);
    vec3 v = vec3(1.0,0.0,0.0);
    float costheta = cos(radians(35.0));
    
    vec3 co = org - center;
    
    float a = dot(dir,v)*dot(dir,v)-costheta*costheta;
    float b = 2.*(dot(dir,v)*dot(co,v) - dot(dir,co)*costheta*costheta);
    float c = dot(co,v)*dot(co,v) - dot(co,co)*costheta*costheta;
    
	// quadratic x^2 = y^2 + z^2
	//float a = dir.y * dir.y + dir.z * dir.z - dir.x * dir.x;
	//float b = org.y * dir.y + org.z * dir.z - org.x * dir.x;
	//float c = org.y * org.y + org.z * org.z - org.x * org.x;
    
	float cap = (s - org.x) / dir.x;
    
	float delta = b * b - 4.0 * a * c;
	if( delta < 0.0 )
		return false;
    
	// 2 roots
	float deltasqrt = sqrt(delta);
	float arcp = 1.0 / a;
	near = (-b - deltasqrt) * 0.5 * arcp;
	far = (-b + deltasqrt) * 0.5 * arcp;
	
	// order roots
	float temp = min(far, near);
	far = max(far, near);
	near = temp;
    
	float xnear = org.x + near * dir.x;
	float xfar = org.x + far * dir.x;
	if( xnear < 0.0 )
	{
		if( xfar < 0.0 || xfar > s )
			return false;
		
		near = far;
		far = cap;
	}
	else if( xnear > s )
	{
		if( xfar < 0.0 || xfar > s )
			return false;
		
		near = cap;
	}
	else if( xfar < 0.0 )
	{
		// The apex is problematic,
		// additional checks needed to
		// get rid of the blinking tip here.
		far = near;
		near = cap;
	}
	else if( xfar > s )
	{
		far = cap;
	}
	    
    //Cone Radius
    float cr = dot(co,co);
    // Inverse Squared Density
    float t1 = near;
    float t2 = far;
    
    float i1 = -(c*t1 + b*t1*t1 + a*t1*t1*t1/3.0);
    float i2 = -(c*t2 + b*t2*t2 + a*t2*t2*t2/3.0);

    density = (i2-i1) / (cr * cr);
    density = min(density, 200.f);
	
    float bValid = 1.0;
    
    density *= bValid;

    float scale =  (3.0 / 4.0) / cr;
    
    density = density*scale;
    
	return true;
}

void transformray (vec3 ro, vec3 rd, mat2 rotationY, vec3 offset, out vec3 outro, out vec3 outrd)
{
	outro = ro + offset;
	outro = vec3(rotationY * outro.xz, outro.y).xzy;
	outrd = vec3(rotationY * rd.xz, rd.y).xzy;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
	// camera
	vec2 q = fragCoord.xy/iResolution.xy;
	vec2 p = -1.0 + 2.0 * q;
	p.x *= iResolution.x/iResolution.y;
	vec3 camro = normalize(vec3(1.0, -0.1, 0.0));
	vec3 w = -camro;
	camro *= 2.5;
	vec3 u = normalize(cross( vec3(0.0, 1.0, 0.0), w ));
	vec3 v = normalize(cross(w,u));
	vec3 camrd = normalize(p.x * u + p.y * v + 1.5 * w);
	fragColor = vec4(0.0);
	
	// rotation
	float angle = 5.0 * iMouse.x / iResolution.x;
	if( iMouse.z < 0.5 )
		angle = iTime + 4.7;
	float ca = cos(angle);
	float sa = sin(angle);
	mat2  m = mat2(ca, -sa, sa, ca);
	
	float far, near;
	vec3 ro, rd;
    float density = 0.0;

	// cone
	transformray(camro, camrd, m, vec3(0, 0.0, 0), ro, rd);
	if (cone (ro, rd, near, far, density))
		//fragColor += vec4(far - max(near, 0.0));
		fragColor += vec4(density);
}

//Book of Shaders Practice and Examples
// MATH: Shaping Functions