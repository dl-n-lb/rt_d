import std.stdio;

// SINGLE FILE WAVEFRONT RAYTRACER
// TODO: IMPLEMENT RAYTRACING OF LAMBERT SPHERES
// TODO: ADD SUPPORT FOR N MATERIALS
// TODO: ADD SUPPORT FOR GENERTIC SHAPES
// TODO: ADD SUPPORT FOR CAMERA AT RUNTIME
// TODO: ACCELERATION STRUCTURE
// TODO: MIS / SHADOW RAYS
// TODO: READ SCENE FROM FILE AT RUNTIME

alias real_t = float;

////////////////////////////////////////////////////////////////////////
///////////////////////////// VECTOR3 //////////////////////////////////
////////////////////////////////////////////////////////////////////////

@nogc struct vec3 {
  real_t x, y, z;

  @nogc pure const vec3 opBinary(string op)(in vec3 rhs) {
    static if (op == "+") {
      return vec3(x + rhs.x, y + rhs.y, z + rhs.z);
    } else if (op == "-") {
      return vec3(x - rhs.x, y - rhs.y, z - rhs.z);
    } else if (op == "*") {
      return vec3(x * rhs.x, y * rhs.y, z * rhs.z);
    } else if (op == "/") {
      return vec3(x / rhs.x, y / rhs.y, z / rhs.z);
    } else {
      assert(0, "Invalid operation");
    }
  }

  @nogc pure const vec3 opBinary(string op)(in real_t rhs) {
    static if (op == "+") {
      return vec3(x + rhs, y + rhs, z + rhs);
    } else if (op == "-") {
      return vec3(x - rhs, y - rhs, z - rhs);
    } else if (op == "*") {
      return vec3(x * rhs, y * rhs, z * rhs);
    } else if (op == "/") {
      return vec3(x / rhs, y / rhs, z / rhs);
    } else {
      assert(0, "Invalid operation");
    }
  }

  @nogc pure const vec3 opBinaryRight(string op)(in real_t lhs) {
    static if (op == "+") {
      return vec3(x + lhs, y + lhs, z + lhs);
    } else if (op == "-") {
      return vec3(lhs - x, lhs - y, lhs - z);
    } else if (op == "*") {
      return vec3(x * lhs, y * lhs, z * lhs);
    } else {
      assert(0, "Invalid operation");
    }
  }
}

@nogc static pure real_t length_squared(in vec3 v) {
  return v.x*v.x + v.y*v.y + v.z*v.z;
}

@nogc static pure real_t length(in vec3 v) {
  import std.math : sqrt;
  return sqrt(v.length_squared);
}

@nogc static pure vec3 norm(in vec3 v) {
  return v / v.length;
}

@nogc static pure vec3 cross(vec3 lhs, vec3 rhs) {
  return vec3(
      lhs.y * rhs.z - lhs.z * rhs.y,
      lhs.z * rhs.x - lhs.x * rhs.z,
      lhs.x * rhs.y - lhs.y * rhs.x
  );
}


////////////////////////////////////////////////////////////////////////
///////////////////////////// RANDOM ///////////////////////////////////
////////////////////////////////////////////////////////////////////////

// FIXME: THE ENEMY OF MANKIND IS GLOBAL STATE
import std.random: Random, uniform01;
auto rng = Random(42);

@nogc static rand_real(in real_t min, in real_t max) {
  return (max - min) * rng.uniform01!real_t() + min;
}

@nogc static vec3 random_vec_in_unit_disk() {
  while(true) {
    vec3 p = vec3(rand_real(-1, 1), rand_real(-1, 1), 0);
    if (p.length >= 1) continue;
    return p;
  }
}

@nogc static pure vec3 color_cast(in vec3 color) {
  import std.math.traits;
  import std.math : sqrt;
  vec3 color_fixed = color;
  if (color.x.isNaN || color.y.isNaN || color.z.isNaN) {
    color_fixed = vec3(1, 0, 1);
  }
  import std.algorithm.comparison : clamp;
  return vec3(
    cast(int)(clamp(sqrt(color_fixed.x), 0, 0.999) * 255),
    cast(int)(clamp(sqrt(color_fixed.y), 0, 0.999) * 255),
    cast(int)(clamp(sqrt(color_fixed.z), 0, 0.999) * 255),
  );
}


////////////////////////////////////////////////////////////////////////
/////////////////////////////// RAY ////////////////////////////////////
////////////////////////////////////////////////////////////////////////

@nogc struct ray_t{
  vec3 o, d;
}

@nogc struct ray_info_t {
  ray_t ray;
  uint x, y;
  vec3 attenuation;
}

alias mat_idx_t = size_t;

@nogc struct hit_info_t {
  ray_t in_ray;
  vec3 p;
  mat_idx_t mat_idx;
}



////////////////////////////////////////////////////////////////////////
///////////////////////////// SHAPES ///////////////////////////////////
////////////////////////////////////////////////////////////////////////

alias object_index = uint; 

// BVH IS CONSTRUCTED OF THESE
struct aabb_t {
  vec3 min, max;
}

// FIXME: how to tell whether its a leaf?
// ig just nullptr for left and right
// should work fine
struct bvh_node_t {
  bvh_node_t* left, right;
  aabb_t box;
  object_index idx;
}

// TODO: HOW TO GENERIC TYPE OBJECT OR SMTH IDK
// IDEA: have an index like 0b1001011101
// where first byte (or smth) represents the type of shape (1 byte gives 256 shapes)
// and remaining bytes give index into a static array (3 bytes gives *enough* of each shape)
// this may not be enough triangles however!!
// this means leaf nodes would be ones where the first byte doesnt represent an aabb ? maybe
struct sphere_t {
  vec3 pos;
  float r;
}


// TODO:
// somehow use mixins to generate all lists of all object types
// array size can be determined at compile since the scene is loaded
// and parsed in an enum
// this can be used to determine N (max of all lengths - may reconsider for triangles)
// actually might be best to do a separate trianglemesh which contains those
// since their numbers will be wayy higher
struct object_list_t(int N, shape_ts...) {
  static foreach(shape_t; shape_ts) {
    mixin(shape_t.stringof ~ "[" ~ N.stringof ~ "] " ~ shape_t.stringof ~ "s;");
  }
}



////////////////////////////////////////////////////////////////////////
/////////////////////////// MATERIALS //////////////////////////////////
////////////////////////////////////////////////////////////////////////

alias material_index = uint;
struct lambert_t {
  vec3 albedo;
}

struct metal_t {
  vec3 albedo;
  float roughness;
}


// TODO: LOAD SCENE
auto load_scene_from_json(JSONValue* json_scene) {
  
}

////////////////////////////////////////////////////////////////////////
///////////////////////////// CAMERA ///////////////////////////////////
////////////////////////////////////////////////////////////////////////

// THE CAMERA HOLDS ALL OF THE INFORMATION TO GENERATE RAYS
// TODO: COMPLETE THIS CLASS
@nogc struct camera_t {
  real_t focus_distance, aperture, fov, samples_per_pixel;

  vec3 lower_left;
  vec3 i, j; // unit vectors across and up
  vec3 lookfrom;

  this(in vec3 lookfrom, in vec3 lookat, in vec3 up, 
       in real_t aspect_ratio, in real_t aperture,
       in real_t fov, in real_t spp) {
    import std.math.constants : PI;
    import std.math : tan;
    real_t half_theta = fov * PI / 360;
    auto vp_h = 2 * tan(half_theta);
    auto vp_w = vp_h * aspect_ratio;

    auto w = (lookfrom - lookat).norm;
    auto u = w.cross(up).norm;
    auto v = w.cross(u);

    this.lookfrom = lookfrom;
    this.focus_distance = (lookfrom-lookat).length;
    this.i = focus_distance * vp_w * u;
    this.j = focus_distance * vp_h * v;
    this.lower_left = lookfrom - i / 2 - j / 2 - focus_distance * w;
    this.aperture = aperture;
    this.fov = fov;
    this.samples_per_pixel = spp;
  }

  private real_t _u = -1;
  private real_t _v = -1;
  private real_t _sample = 0;

  // TODO: GENERATE THE NEXT RAY TO BE CAST - NON PARALLEL IS FINE
    @nogc ray_info_t generate_ray(int width, int height) {
    if (_u == -1) _u = -1.0 / real_t(width);
    if (_v == -1) _v = 0;
    _u += 1.0 / real_t(width);
    if (_u >= 1.0 - 1.0/real_t(width)) {
      _u = 0;
      _v += 1.0/ real_t(height);
      if (_v >= 1.0 - 1.0/real_t(height)) {
        _u = -1;
        _v = -1;
        ++_sample;
      }
    }
    vec3 rd = aperture / 2 * random_vec_in_unit_disk();
    vec3 offset = i * rd.x + j * rd.y;
    return ray_info_t(
      ray_t(
        lookfrom + offset, 
        lower_left + _u * i + _v * j - lookfrom - offset
      ), 
      cast(uint)(_u * width), cast(uint)(_v * height),
      vec3(1, 1, 1)
    );   
  }

  @nogc bool can_generate_rays() {
    
    return (_u < 0.999999 || _v < 0.999999) && samples_per_pixel > _sample;
  }
}

import std.json: JSONValue, JSONType;

auto from_json(JSONValue cam_json, uint width, uint height) {
  import core.stdc.stdlib;
  vec3 lookfrom, lookat = vec3(0, 0, 0), up = vec3(0, 1, 0);
  real_t aperture = 0.0, fov = 45, aspect = real_t(width)/height, spp = 1;
  if ("Camera" !in cam_json) { writeln("User must provide camera settings!"); exit(1); }
  if (const(JSONValue)* json_lookfrom = "Position" in cam_json["Camera"]) {
    auto lf = json_lookfrom.array;
    lookfrom = vec3(lf[0].get!float, lf[1].get!float, lf[2].get!float);
  } else {
    writeln("Camera file must provide a position!"); exit(1);
  }
  if (const(JSONValue)* json_lookat = "LookAt" in cam_json["Camera"]) {
    auto la = json_lookat.array;
    lookat = vec3(la[0].get!float, la[1].get!float, la[2].get!float);
  }
  if (const(JSONValue)* json_up = "Up" in cam_json["Camera"]) {
    auto u = json_up.array;
    lookat = vec3(u[0].get!float, u[1].get!float, u[2].get!float);
  }
  if (const(JSONValue)* j = "Aperture" in cam_json["Camera"]) {
    auto json_aperture = *j;
    aperture = json_aperture.get!float;
  }
  if (const(JSONValue)* j = "fov" in cam_json["Camera"]) {
    auto json_fov = *j;
    fov = json_fov.get!float;
  }
  if ("Settings" !in cam_json) { writeln("User must provide render settings!"); exit(1); }
  if (const(JSONValue)* j = "Samples" in cam_json["Settings"]) {
    spp = (*j).get!uint;
  }
  return camera_t(lookfrom, lookat, up, aspect, aperture, fov, spp);
}

enum output_type {
  PPM,
  PNG
}

// FIXME: HOW TO GET PIXEL DATA FROM THIS SETUP
// CAN ADD A size_t TO ray_info_t AND hit_info_t BUT THIS IS POLLUTING THE CACHE
// IS THERE A BETTER WAY PERHAPS A SEPARATE ARRAY THAT STORES THEM THATS UPDATED EACH LOOP
// THE OBFUSCATION OF THE CODE WORTH A SIZE REDUCTION OF 4-8 BYTES?
struct app_t(int width, int height, uint batch_size) {
  camera_t camera;

  ray_info_t[batch_size] ray_queue;
  hit_info_t[batch_size] hit_queue;
  size_t[batch_size] curr_ray_pixels;

  vec3[width*height] output;

  this(in camera_t camera) {
    this.camera = camera;
  }

  // TODO: THIS IS THE SCENE INTERSECTION
  // SINCE THIS IS MOSTLY BVH TRAVERSAL IT SHOULD BE SIMPLE TO PARALLELIZE
  // HOPEFULLY
  private @nogc void cast_rays() {

  }

  // TODO: THIS IS THE MATERIAL CALCULATION
  // SOME SORT OF META-CODE TO GENERATE EACH MATERIAL AS A PARALLEL LOOP MAYBE
  // BUT THIS MEANS THAT MOST OF THE DATA ISN'T OPERATED ON AT ONCE
  // PERHAPS ONE MATERIAL PER THREAD
  // FOR CONCURRENCY OF SORTS
  private @nogc void generate_hit_rays() {

  }

  // TODO: THIS SECTION IS THE PART WHICH IS HARD TO PARALLELIZE
  private @nogc size_t flatten_ray_queue() {
    return batch_size;
  }

  private @nogc void refill_ray_queue(size_t curr_end_index) {

  }

  @nogc void draw(uint max_depth) {
    // TODO: generate first queue of rays
    // FIXME THIS LOOP IS DISABLED FOR NOW
    for (int i = 0; i < batch_size; ++i) {
      ray_queue[i] = camera.generate_ray(width, height);
    }

    while (camera.can_generate_rays()) {
      cast_rays();
      generate_hit_rays();
      size_t end_index = flatten_ray_queue();
      refill_ray_queue(end_index);
      camera.generate_ray(width, height);
    }
    // temp
    for (int i = 0; i < width*height; ++i) {
      vec3 rd = camera.generate_ray(width, height).ray.d.norm;
      real_t t = 0.5 * (rd.y + 1);
      output[i] = ((1 - t) * vec3(1, 1, 1) + t * vec3(0.4, 0.6, 1.0)).color_cast;
    }
  }

  void write_file(string filename, output_type type) {
    switch(type) {
      case output_type.PPM:
        auto f = File(filename, "w");
        f.writeln("P3");
        f.writeln(width, " ", height, "\n255");
        foreach(vec3 pixel; output) {
          f.writeln(pixel.x, " ", pixel.y, " ", pixel.z);
        }
        break;
      case output_type.PNG:
        import dlib.image;
        SuperImage img = image(width, height);
        for (int i = 0; i < width*height; i++) {
          vec3 p = output[i];
          img[i % width, i / width] = Color4f(p.x, p.y, p.z);
        }
        img.savePNG(filename);
        break;
      default: break;
    }
  }
}

import std.file;
import std.json;

void main()
{
  enum scenefile = import("source/scene.json").parseJSON();

  // TODO: camera generated from json later
  // TODO: Get json from exe arguments later
  //auto camera_file = readText("camera.json").parseJSON();
  auto camera_file = import("camera.json").parseJSON();
  auto cam = from_json(camera_file, 800, 600);
  uint max_depth = 3;
  if (const(JSONValue)* j = "MaxDepth" in camera_file["Settings"]) {
    max_depth = (*j).get!uint;
  }

  auto app = app_t!(800, 600, 500)(cam);
  app.draw(max_depth);
  app.write_file("image", output_type.PPM);
}
