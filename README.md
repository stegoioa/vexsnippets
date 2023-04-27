<pre>

# vexsnippets
Useful VEX Code Snippets

VEX Code Snippets
This repository contains a collection of VEX code snippets that are commonly used in Houdini workflows. The purpose of this repository is to provide a quick and easy way for Houdini users to integrate these code snippets into their VEX pipelines.

### rand_removepoint1

// removes random point with certain probability 
// run over points, create spare parameters
removepoint(0, rand(@ptnum + ch("seed")) < ch("probability")?@ptnum:-1);
 

rand_removepoint2

// removes random point with certain probability 
// run over points, create spare parameters
float prob = chf('probability');
float seed = chf('seed');
float u = rand(set(@elemnum%666, @elemnum/666, seed));
if (prob > u){
removepoint(geoself(), @ptnum);
}
 

rand_removeprim1

// removes random primitive with certain probability 
// run over primitives, create spare parameters
removeprim(0, rand(@primnum + ch("seed")) < ch("probability")?@primnum:-1, 1);
 

rand_removeprim2

// removes random primitive with certain probability 
// run over primitives, create spare parameters
float prob = chf('probability');
float seed = chf('seed');
float u = rand(set(@elemnum%666, @elemnum/666, seed));
if (prob > u){
removeprim(geoself(), @primnum, 1);
}
 

rand_pscale1

// create spare parameters
float rand = rand(@ptnum+ch('seed'));
float pow = ch('power'); // default is 1
@pscale = fit01(pow(rand, pow), ch('min'), ch('max'));
 

rand_color_ramp1
Source

v@Cd = vector(chramp('color', rand(@ptnum+ch('seed'))));
 

rand_rotation_axis1

// set your axis or use existing from normals
// create spare parameters
float minmaxangle = chf('minmax_angle')*PI*2;
float angle = fit01(rand(@ptnum+ch('seed')), -minmaxangle, minmaxangle);
@rot = quaternion(angle, chv('axis'));
int usenormals = chi('use_normals');
if (usenormals == 1){ // 0-1 switch
    @rot = quaternion(angle, v@N);
}
 

rand_orient1
Source

// create spare parameters
float mm_x = ch('minmax_x')*180;
float mm_y = ch('minmax_y')*180;
float mm_z = ch('minmax_z')*180;
float angl_x = radians(fit01(rand(@ptnum+ch('seed')), -mm_x, mm_x));
float angl_y = radians(fit01(rand(@ptnum+ch('seed')), -mm_y, mm_y));
float angl_z = radians(fit01(rand(@ptnum+ch('seed')), -mm_z, mm_z));
vector angle = set(angl_x, angl_y, angl_z);
@orient = eulertoquaternion(angle, 0);
 

set_orient1
Source

// create spare parameters
@orient = quaternion(maketransform(normalize(-@P),{0,1,0}));
vector4 pitch = quaternion({1,0,0}*ch('pitch')*PI*2);
vector4 yaw   = quaternion({0,1,0}*ch('yaw')*PI*2);
vector4 roll  = quaternion({0,0,1}*ch('roll')*PI*2);
@orient = qmultiply(@orient, pitch); //change order
@orient = qmultiply(@orient, yaw);
@orient = qmultiply(@orient, roll);
 

vector_dir1

// create spare parameters
@`chs('vec_attrib')` = chv('dir');
 

peak_vex1

// run over points, create spare parameters
float min = ch('min');
float max = ch('max');
float seed = ch('seed');
float randpeak = fit01(rand(@ptnum + seed), min, max);
@P += @N * randpeak;
 

lookat1

// points are looking to the second input
// run over points
@N = normalize(point(1, 'P', 0) - @P);
@up = {0, 1, 0};
 

xz_jitter1
Source

// create spare parameters
@P += curlnoise(@P)*{1,0,1}*fit01(rand(chi('seed')),ch('min'),ch('max'));
 

simple_spherify1
Source

// create spare parameters
@P = lerp(@P, normalize(@P - getbbox_center(0)) * ch('rad') + getbbox_center(0), ch('amt'));
 

noise_gen1

// create spare parameters and check defaults
// can be used for Cd, P noise or mountains
// x- s- o- a- curl- noise functions
float scale = chf('scale'); // default is 1
vector freq = ch('freq');
vector offset = chv('offset');
int turb = chi('turb');
float rough = ch('rough');
float atten = ch('atten'); // default is 1
vector noise = onoise(@P*scale+freq+offset, turb, rough, atten);
noise *= ch('amp');
int invert = chi('invert'); // switch 0-1
if (invert == 1){
    noise *= -1;
    }
int noisetype = chi('cd____p____mountains'); // switch 0-2
if (noisetype == 0){ // Cd noise
    @Cd = fit(noise[chi('color_channel')], -0.5, 0.5, 0, 1); // 0-2
    }
if (noisetype == 1){ // P noise
    @P += noise;
    }
if (noisetype == 2){ // mountains on the ZX plane
    @P.y += abs(noise[1]);
    }
 

noise_normals1

// create spare parameters and check defaults
// can be used for N noise
// x- s- o- a- curl- noise functions
float scale = chf('scale'); // default is 1
vector freq = ch('freq');
vector offset = chv('offset');
int turb = chi('turb');
float rough = ch('rough');
float atten = ch('atten'); // default is 1
vector noise = onoise(@P*scale+freq+offset, turb, rough, atten);
noise *= ch('amp');
@N = normalize(@N + noise);
 

cross_product1

// run over points, create spare parameters
@N = cross(@P, chv('axis'));
@v = -@N;
 

dot_product_group1
Source

// run over primitives, create spare parameters
// use dotgroup1 in bindings > output selection group
float angle = acos(dot(normalize(chv('up_vector')), normalize(@N)));
angle /= (PI / 180);
if (angle < ch('threshold')){
    @group_dotgroup1 = 1;
    }
 

dot_product_density1

// run over points, create spare parameters
// use attribute blur afterward
float angle = acos(dot(normalize(chv('up_vector')), normalize(@N)));
angle /= (PI / 180);
if (angle < ch('threshold')){
    @density = 1; } else { @density = 0;
    }
if (chi('visualize') == 1){ // switch 0-1
    @Cd = @density;
    }
 

zero_out_y1

// zero out points on Y axis
@P = set(@P.x, ch('level'), @P.z);
 

spiralize1

// requires dots or line as an input
// create spare parameters
float freq = ch('freq');
float rad = ch('rad');
@P.x = cos(@ptnum * freq)*@P.y*rad;
@P.z = sin(@ptnum * freq)*@P.y*rad;
@P.y *= ch('height');
 

color_harm_gen1
Source

// create spare parameters
@Cd = rgbtohsv(set(lerp(0.3, 1, point(0, 'pscale', @ptnum)*7), random(@ptnum)*0.7, random(@ptnum+30)*0.3)*1.2);
@Cd = hsvtorgb(@Cd+set(ch('hue'), 0, 0));
 

raycasting1
Source

// intersection between input 0 and 1
// create spare parameters
vector hit_dir = chv('hit_direction');
vector hit_pt;
vector uvw;
int prim_num = intersect(1, @P, hit_dir, hit_pt, uvw);
if (prim_num >= 0){
    @P = hit_pt;
    @N = prim(1, 'N', prim_num);
    @group_rayHitGroup = 1;
    }
 

color_ramp_bbox1
Source

// create spare parameters
vector bb = relbbox(0, @P);
@Cd = vector(chramp('ramp', bb.y));
 

color_correction1
Source

// create spare parameters
v@Cd=hsvtorgb(rgbtohsv(v@Cd)+set(ch('hue'),ch('sat'),ch('val')));
 

midpt_group1
Source

// creates a group with the middle points
i@group_midpt1 = neighbourcount(0, @ptnum) == 2;
 

endpt_group1
Source

// creates a group with the first and the last points
if (@ptnum == 0 || @ptnum == @numpt-1) i@group_pin = 1;
 

endpt_group2
Source

// creates a group with the first and the last points
int first_vtx = primvertex(0, @primnum, 0);
int first_pt = vertexpoint(0, first_vtx);
int last_vtx = primvertex(0, @primnum, primvertexcount(0, @primnum)-1);
int last_pt = vertexpoint(0, last_vtx);
setpointgroup(0, 'endpt', first_pt, 1);
setpointgroup(0, 'endpt', last_pt, 1);
 

second_input_attribute1
Source

// examples
i@number_of_points_input_1 = npoints(0);
v@Cd = point(1, 'Cd', @ptnum);
v@pos_input_1 = point(1, 'P', @ptnum);
v@pos_input_1_1 = v@opinput1_P;
 

normal_orient1

@N = set(0, -1, 0);
 

group_by_attribute1

if(@testattrib==0.666){
@group_testgroup=1;
}
 

number_round1

// create spare parameters
float step = chf("step"); // put 0.01
@testattrib = rint(@pscale/step) * step;
 

quantize_pos1
Source

float ps = chf("scale");
@P.x = rint(@P.x/ps) * ps;
@P.y = rint(@P.y/ps) * ps;
@P.z = rint(@P.z/ps) * ps;
 

median_pt_num1
Source

int npt = npoints(0);
if(npt%2==0)
    i@midpt = (npt/2) -1;
else
    i@midpt = (int(ceil(float(npt)/2))) -1;;
 

points_box1

// creates box with points using vex
// run over numbers
// create spare parameters
addpoint(0, fit01(rand(@elemnum), chv('min'), chv('max')));
 

hex_grid1
Source


create_hex_corners1

// run over details

vector  Global_Start_Pt, Hexagon_Start_Pt;
vector  Apex_Base_Plot, Apex_Current_Plot, Current_Hex_Center;
float   Hexagon_Radius;
int     Hexagons_Per_Row, Number_Of_Rows;
float   Angle;
vector  Axis;
matrix3 MTX_For_Rotation;
matrix  MTX_For_Translation;
int     Count, Count_A, Count_B;
int     Pt_For_Set;

//*******************************************************

Global_Start_Pt      = set(0, 0, 0);
Hexagon_Start_Pt     = chv('starting_point');
Hexagon_Radius       = chf('radius');
Hexagons_Per_Row     = chi('elements_per_row');
Number_Of_Rows       = chi('number_of_rows');

i@Total_Hexagons     = Hexagons_Per_Row * Number_Of_Rows;

//*******************************************************

float First_Adj, Second_Adj, First_Hyp;

First_Hyp  = Hexagon_Radius;
First_Adj  = cos(radians(30)) * First_Hyp;
Second_Adj = First_Adj / tan(radians(30));

//*******************************************************

Axis               = set(0,1,0);
Current_Hex_Center = Global_Start_Pt;

for(Count_B = 0; Count_B < Number_Of_Rows; Count_B++)
    {

    Current_Hex_Center = Global_Start_Pt;
    Current_Hex_Center.z +=  ((First_Adj * 2) * Count_B) - First_Adj;

    for(Count_A = 0; Count_A < Hexagons_Per_Row; Count_A++)
        {

        Apex_Base_Plot   = Global_Start_Pt;
        Apex_Base_Plot.x += Hexagon_Radius;

        Current_Hex_Center.x = Global_Start_Pt.x + (Second_Adj * Count_A);

        if( (Count_A % 2) == 0) { Current_Hex_Center.z += First_Adj; }           
        if( (Count_A % 2) == 1) { Current_Hex_Center.z -= First_Adj; }          

        for(Count = 0; Count < 6; Count++)
            {
            MTX_For_Rotation    = ident();
            MTX_For_Translation = ident();

            Angle = radians(60 * Count);

            rotate(MTX_For_Rotation, Angle, Axis);
            translate(MTX_For_Translation, (Current_Hex_Center + Hexagon_Start_Pt));

            Apex_Current_Plot = Apex_Base_Plot;
            Apex_Current_Plot *= MTX_For_Rotation;
            Apex_Current_Plot *= MTX_For_Translation;

            Pt_For_Set = addpoint(geoself(), Apex_Current_Plot);
            setpointattrib(geoself(), 'Hex_Num', Pt_For_Set, ((Count_B * Hexagons_Per_Row) + Count_A), 'set');
            }

        }  
    }

connect_to_polylines_or_polygons1

// run over details

int Count, Count_A;
int Total, Prim_Num;
int Polytype;

Polytype = chi('polytype'); // 0-1 switch
Total = detail(geoself(), 'Total_Hexagons');

if (Polytype == 0)
    {
    for(Count = 0; Count < Total; Count++)
        {
        Prim_Num = addprim(geoself(), 'polyline');         
        addvertex(geoself(), Prim_Num, (Count * 6)    );
        addvertex(geoself(), Prim_Num, (Count * 6) + 1);
        addvertex(geoself(), Prim_Num, (Count * 6) + 2);
        addvertex(geoself(), Prim_Num, (Count * 6) + 3);
        addvertex(geoself(), Prim_Num, (Count * 6) + 4);
        addvertex(geoself(), Prim_Num, (Count * 6) + 5); 
        addvertex(geoself(), Prim_Num, (Count * 6)    );            
        }
    }

if (Polytype == 1)
    {
    for(Count = 0; Count < Total; Count++)
        {
        Prim_Num = addprim(geoself(), 'poly');         
        addvertex(geoself(), Prim_Num, (Count * 6)    );
        addvertex(geoself(), Prim_Num, (Count * 6) + 1);
        addvertex(geoself(), Prim_Num, (Count * 6) + 2);
        addvertex(geoself(), Prim_Num, (Count * 6) + 3);
        addvertex(geoself(), Prim_Num, (Count * 6) + 4);
        addvertex(geoself(), Prim_Num, (Count * 6) + 5);            
        }
    }
 

hex_grid2

// hexagonal point grid generator
// run over details
int geo = geoself();
int size_y = chi('size_y');
int size_x = chi('size_x');
int inv_order = chi('inv_order');
for (int j = 0; j < size_y; j++){
    float shift_x = ((j + inv_order) % 2) * 0.5;
    float shift_z = j * (sqrt(3) * 1) * 0.5;
    for (int i = 0; i < size_x; i++){
        vector pos = set(i + shift_x, 0, shift_z);
        addpoint(-1, pos);
    }
}
 

bend_curve1
Source

float amount = ch("amount");
float gradient = (float)@ptnum/(float)(@numpt-1);
@P.y += pow((2 * gradient) -1, 2) * amount;
@P.y -= amount;
 

bend_curve2
Source

// create spare parameters
// input is a line
// resample SOP on the line for more points
// activate curveu attrib in resample SOP
// use transform SOP with Y translate -$YMAX for 0 yposition
@curveu=chramp("ramp",@curveu);
float bendamt = chf("bend_amount");
vector benddir = chv("bend_direction");
@P+= benddir * bendamt * @curveu;
 

Parabola formula to create wires
Extra info | Source


parabolic_tension1

// minimum 3 points is required
// create spare parameters
float ptn = @ptnum;
float npt = @numpt;
float curveu = @ptnum / (npt-1);
float paraby = 1 - (pow((2 * curveu-1), 2));
@P.y -= paraby * ch("tension");
 

color_slicer1
Source

// create spare parameters
// use color ramp
// remove_even_odd, even_odd, show_color are 0-1 toggles

float cursegmentmax, cursegmentmin, evenodd, grey;
int pointsinslice[], pt;

int slices = chi('slices');
float reo = chi('remove_even_odd');
float eo = chi('even_odd');
int seed = chi('seed');

vector bb = getpointbbox_size(0);
vector bbmin = getpointbbox_min(0);
float segment = bb.y / slices;

for (int i = 1; i <= slices; i++){
    grey = fit(i, 0, slices, 0, rand(i*seed));
    
    cursegmentmax = (segment * i) + bbmin.y;
    cursegmentmin = cursegmentmax - segment;
    
    if(@P.y >= cursegmentmin && @P.y <= cursegmentmax){
        @Cd = set(grey, grey, grey);
        i@idgrp = i;
        push(pointsinslice, @ptnum);
    }
    
    if(reo == 1){
        evenodd = @idgrp%2;
        
        if(eo == 0){
            if(evenodd == 0){
                foreach(pt; pointsinslice){
                    removepoint(0, pt);
                }
            }
        }else{
            if(evenodd == 1){
                foreach(pt; pointsinslice){
                    removepoint(0, pt);
                    }
                }
            }
        }
    }

float showcolor = chi('show_color');
if (showcolor == 1){
    @Cd = chramp('color_ramp', @Cd.r);
}
 

Replace a prim with a point
Source


prim_to_point1

// run over primitives
addpoint(0, @P);
removeprim(0, @primnum, 1);

prim_to_point_n1

// with point normals
// run over primitives
int newpoint = addpoint(0, @P);
setpointattrib(0, “N”, newpoint, @N);
@up = set(0, 0, 1); 
removeprim(0, @primnum, 1);
 

object_center1
Source


// centroid computed from the second input
// run over detail
vector centro = getbbox_center(1);
addpoint(0, centro);
 

String to number convertion
Source


str_to_num_conv1

// simple string to number convertion
f@out=atof("12.5");
i@outInt=atoi("12");

str_to_num_conv2

// string to number convertion with other symbols
s@in="frgh1234.2akjf";
f@out=atof(re_replace(r"[^.0-9]+","",@in));

str_to_num_conv3

// parsing string values with multiple numbers
s@in="123hjh765.4hsdjfh12hh";
string array[]=re_split(r"[^.0-9]+",@in);
 

Procedural subdivision curves
Source


procedural_subd_curves_x1

// connect equal amount of points to inputs 0 and 1
// works best along X-axis
// create spare parameters
// add Resample SOP after
float o = chf("Offset");
float rando = fit01(rand(@ptnum + chi("Seed")), -o, o) * chf("Random_Offset_Amp");
o += rando;
float s = chf("Curve_Softness");
vector npos = point(1, "P", @ptnum);
vector pos = v@P;
vector d = npos - pos;
vector p2 = pos + set(clamp((o - s*o), 0, 1) * d.x, 0, 0);
vector p3 = pos + set(o * d.x, 0, d.z * (1.0 - s));
vector p4 = pos + set(clamp((o + s*o), 0, 1) * d.x, d.y, d.z);
int npt2 = addpoint(0, p2);
int npt3 = addpoint(0, p3);
int npt4 = addpoint(0, p4);
int npt5 = addpoint(0, npos);
int newprim = addprim(0, "polyline", @ptnum, npt2, npt3, npt4);
addvertex(0, newprim, npt5);

procedural_subd_curves_y1

// connect equal amount of points to inputs 0 and 1
// works best along Y-axis
// create spare parameters
// add Resample SOP after
float o = chf("Offset");
float rando = fit01(rand(@ptnum + chi("Seed")), -o, o) * chf("Random_Offset_Amp");
o += rando;
float s = chf("Curve_Softness");
vector npos = point(1, "P", @ptnum);
vector pos = v@P;
vector d = npos - pos;
vector p2 = pos + set(0, clamp((o - s*o), 0, 1) * d.y, 0);
vector p3 = pos + set(d.x * (1.0 - s), o * d.y , 0);
vector p4 = pos + set(d.x, clamp((o - s*o), 0, 1) * d.y, d.z);
int npt2 = addpoint(0, p2);
int npt3 = addpoint(0, p3);
int npt4 = addpoint(0, p4);
int npt5 = addpoint(0, npos);
int newprim = addprim(0, "polyline", @ptnum, npt2, npt3, npt4);
addvertex(0, newprim, npt5);
 

Connects any found points; supports second input detection
Source


rand_nearby_connect1

// create spare parameters
int read=0;
if(npoints(1)>0){
read=1;
}
int pts[]=pcfind(read,"P",v@P,chf("maxLineLength"),chi("maxFindCount"));
int l=len(pts);
float fl=float(l);
int randomConnect=chi("randomConnect");
int rander=pts[int(random(@ptnum*fl+randomConnect)*fl*fl) % l];
vector curPos=attrib(read,"point", "P", rander);
int to=addpoint(0,curPos);
addprim(0,"polyline",@ptnum,to);

line_generator1

// create spare parameters
int drawLine(vector f; vector t){
    int base=addpoint(0,f);
    int to=addpoint(0,t);
    int line=addprim(0,"polyline",base,to);
    return line;
}
int read=0;
if(npoints(1)>0){
read=1;
}
int maxLines=chi("maxLineCount");
float minLen=chf("minLineLength");
int pts[]=pcfind(read,"P",v@P,chf("maxLineLength"),chi("maxFindCount"));
int randomConnect=chi("randomConnect");
int keepPointCount=min(1, max(0,chi("keepPointCount")));
int runner=0;
vector curPos;
int pt;
if(randomConnect == 0){
    for(int x=0; x<len(pts);++x){
        pt=pts[x];
        if(runner > maxLines){ break; }
        curPos=attrib(read,"point", "P", pt);
        if(length(curPos-v@P)>minLen && (@ptnum<pt || read)){
            if(keepPointCount){
                int to=pt;
                if(read){
                    to=addpoint(0,curPos);
                }
                addprim(0,"polyline",@ptnum,to);
            }else{
                drawLine(v@P,curPos);
            }
            runner++;
        }
    }
}else{
    int l=len(pts);
    float fl=float(l);
    int rander=pts[int(random(@ptnum*fl+randomConnect)*fl*fl) % l];
    curPos=attrib(read,"point", "P", rander);
    if(keepPointCount){
        int to=rander;
        if(read){
            to=addpoint(0,curPos);
        }
        addprim(0,"polyline",@ptnum,to);
    }else{
        drawLine(v@P,curPos);
    }
}
 

Circle Pattern
Source


circle_generator1

// run over detail
// create spare parameters
int sample = chi('sample');
float radius = ch('radius');
vector center = chv('center');
// circle formula
float theta = 0;
float step = PI*2 / (float)sample;
float x,z;
vector pos;
while(theta < PI*2)
{
    x = center.x + cos(theta)*radius;
    z = center.z + sin(theta)*radius;
    pos = set(x, center.y, z);
    addpoint(0, pos);
    theta += step;
}

line_generator1

// create spare parameters
int line = addprim(0, 'polyline');
int const_offset = chi('offset');
float phase = ch('phase');
// additional offset
float fit_ptnum = fit(@ptnum, 0, npoints(0), 0, PI*phase);
int sin_offset = (int)(sin(fit_ptnum) * ch('sin_offset'));
// find connection neighbour
int neigh = (@ptnum + const_offset + sin_offset)%npoints(0);
// create the connection line
addvertex(0, line, @ptnum);
addvertex(0, line, neigh);

color_ramp1

int pts[] = primpoints(0, @primnum);
float fit; vector color;
for (int i=0; i < len(pts); i++)
{
	fit = fit(i, 0, len(pts), 0, 1);
	fit = chramp('color', fit);
	color = set(fit, fit, fit);
	setpointattrib(0, "Cd", pts[i], color, "set");
}
 

Attribute name as a string from channel
Source


attrib_name1

// create spare parameters
string name = chs("attributeName");
int newatt = addattrib(0, "point", name, "...");
 

range_filter1

// create spare parameters
i@everyother = set((@ptnum+chi('offset'))%chi('delete_every')==chi('inv_order'));
if(@everyother){
    removepoint(geoself(), @ptnum);
    }
 

curveu_vex1
Source

// run over points
int prvtx = vertexprimindex(0, @vtxnum);
@curveu = (float)prvtx / (@numvtx - 1);
 

minpos_lerp1

// blends between current and minimum position
// run over points
@P = lerp(@P, minpos(1, @P), ch('amount'));
 

procedural_line1
Source

// create spare parameters
// run over detail
float length = ch('length');
int pointcount = chi('points');
vector dir = chv('direction');
vector startpos = chv('start_point_position');
dir = normalize(dir);
int points[];
resize(points, pointcount);
f@stepvalue = length/(float)(pointcount-1);
for (int i = 0; i < pointcount; i++)
{
    vector pointpos = dir * (@stepvalue*i) + startpos;
    int currentID = addpoint(geoself(), pointpos);
    points[i] = currentID;
}
addprim(geoself(), "polyline", points);
 

recursive_subdivision1
Source

// create spare parameters
// run over primitives
float x = chf("x_amount");
float y = chf("y_amount");
int pnts[] = primpoints(0, @primnum);
vector pos_0, pos_1, pos_2, pos_3;
vector pos_a, pos_b, pos_c, pos_d, pos_e, pos_f;
pos_0 = point(0, "P", pnts[0]);
pos_1 = point(0, "P", pnts[1]);
pos_2 = point(0, "P", pnts[2]);
pos_3 = point(0, "P", pnts[3]);
pos_a = lerp(pos_0, pos_1, y);
pos_b = lerp(pos_1, pos_2, x);
pos_c = lerp(pos_3, pos_2, 1-y);
pos_d = lerp(pos_0, pos_3, x);
pos_e = lerp(pos_d, pos_b, y);
pos_f = lerp(pos_d, pos_b, 1-y);
int pnt_a = addpoint(0, pos_a);
int pnt_b = addpoint(0, pos_b);
int pnt_c = addpoint(0, pos_c);
int pnt_d = addpoint(0, pos_d);
int pnt_e = addpoint(0, pos_e);
int pnt_f = addpoint(0, pos_f);
addprim(0, "poly", pnts[0], pnt_a, pnt_e, pnt_d);
addprim(0, "poly", pnt_a, pnts[1], pnt_b, pnt_e);
addprim(0, "poly", pnt_b, pnts[2], pnt_c, pnt_f);
addprim(0, "poly", pnt_d, pnt_f, pnt_c, pnts[3]);
removeprim(0, @primnum, 0);
 

vein_noise1
Source

vector npos = v@P/1. + set(0., 666., 0.);   // Noise input 3D position
float namp = 1.;                            // namp (Noise amplitude)
float nval = 0., nweight = 0.;              // Init nval (Noise output value), and nweight (Used to normalize octaves)
int oct = 9;                                // Number of Octaves
for( int i = 0; i < oct; i++ )    {
    float __nval = fit(abs(-0.5+noise(set(npos.x,npos.y,npos.z,f@Time))), 0.0, 0.1, 1., 0.);
    nval += __nval * namp;                  // Amplitude
    nweight += namp;                        // Accumulate weight
    npos *= 2.132433;                       // Lacunarity
    namp *= 0.666;                          // Roughness
}
v@Cd = 1 - pow(nval / nweight, 0.8765);     // Visualize Noise Output
 

remove_overlap_prim1
Source

// removes overlapping primitives
int prim_points[] = primpoints(geoself(), @primnum );
vector pos_accum = 0;

for (int i = 0; i < len(prim_points); i++){
    pos_accum += attrib(0, 'point', 'P', prim_points[i]);
}

pos_accum /= len(prim_points);

int xyz_prim;
vector xyz_uv;

float xyz_dist = xyzdist(0, pos_accum, xyz_prim, xyz_uv);

if (xyz_prim > @primnum && xyz_dist < 0.001)
    removeprim(geoself(), @primnum, 1);
else if ( xyz_prim < @primnum && xyz_dist < 0.001)
    removeprim(geoself(), xyz_prim, 1);
    
    <pre>
