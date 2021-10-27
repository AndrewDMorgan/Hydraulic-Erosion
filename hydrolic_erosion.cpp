#include "image.h"  // contains other includes

#define scale 150  // i belive this is the height scaler for the world which determins the rate at which the dropplets accelerate down hill


// a particle class for water droplets
class Particle
{
    // the particle class has no private variables or methods
    public:

        float2 pos;
        float2 speed;

        float volume;
        float sediment;

        // setting both speed and position
        Particle (float2 pos_, float2 speed_)
        {
            speed = speed_;
            pos = pos_;

            sediment = 0;
            volume = 1;
        }

        // for if you dont set speed
        Particle (float2 pos_)
        {
            speed = float2(0);
            pos = pos_;

            sediment = 0;
            volume = 1;
        }
};


// gets the surface normal of the height map at a given position
float3 SurfaceNormal(int i, int j, float2 map_size, Array <float>& height_map)
{
    float3 n = float3(0.);
    n = float3(0.15, 0.15, 0.15) * normalize(float3(scale*(height_map(i,j)-height_map(min(i + 1, (int) map_size.x), j)), 1.0, 0.0));

    n += float3(0.15, 0.15, 0.15) * normalize(float3(scale*(height_map(max(i - 1,0),j)-height_map(i,j)), 1.0, 0.0));
    n += float3(0.15, 0.15, 0.15) * normalize(float3(0.0, 1.0, scale*(height_map(i,j)-height_map(i,min(j + 1, (int) map_size.y)))));
    n += float3(0.15, 0.15, 0.15) * normalize(float3(0.0, 1.0, scale*(height_map(i,max(j - 1,0))-height_map(i,j))));

    n += float3(0.1, 0.1, 0.1) * normalize(float3(scale*(height_map(i,j)-height_map(min(i + 1, (int) map_size.x),min(j + 1, (int) map_size.y)))/sqrt(2), sqrt(2), scale*(height_map(i,j)-height_map(min(i + 1, (int) map_size.x),min(j + 1, (int) map_size.y)))/sqrt(2)));
    n += float3(0.1, 0.1, 0.1) * normalize(float3(scale*(height_map(i,j)-height_map(min(i + 1, (int) map_size.x),max(j - 1,0)))/sqrt(2), sqrt(2), scale*(height_map(i,j)-height_map(min(i + 1, (int) map_size.x),max(j - 1,0)))/sqrt(2)));
    n += float3(0.1, 0.1, 0.1) * normalize(float3(scale*(height_map(i,j)-height_map(max(i - 1,0),min(j + 1, (int) map_size.y)))/sqrt(2), sqrt(2), scale*(height_map(i,j)-height_map(max(i - 1,0),min(j + 1, (int) map_size.y)))/sqrt(2)));
    n += float3(0.1, 0.1, 0.1) * normalize(float3(scale*(height_map(i,j)-height_map(max(i - 1,0),max(j - 1,0)))/sqrt(2), sqrt(2), scale*(height_map(i,j)-height_map(max(i - 1,0),max(j - 1,0)))/sqrt(2)));

    return n;
}


// 1d perlin noise (aka just splines)
float perlin1D(Array <float>& noise, int index_x, float h, float x)
{
    float nx = x / h;
    int p1 = (int) nx;
    int p2 = p1 + 1;
    int p3 = p2 + 1;
    int p0 = p1 - 1;
    
    float t = nx - p1;
    float tt = t * t;
    float ttt = tt * t;
    float ttt3 = ttt * 3;
    
    float q1 = -ttt + 2*tt - t;
    float q2 = ttt3 - 5*tt + 2;
    float q3 = -ttt3 + 4*tt + t;
    float q4 = ttt - tt;

    float ty = 0.5 * (noise(index_x, p0) * q1 + noise(index_x, p1) * q2 + noise(index_x, p2) * q3 + noise(index_x, p3) * q4);

    return ty;
}

// splining the splines to get 2d splines
float perlin2D(Array <float>& noise, float h, float x, float y)
{
    float nx = x / h;
    int i = (int) nx;

    float ty1 = perlin1D(noise, i - 1, h, y);
    float ty2 = perlin1D(noise, i    , h, y);
    float ty3 = perlin1D(noise, i + 1, h, y);
    float ty4 = perlin1D(noise, i + 2, h, y);

    Array <float> n_array = Array <float> (1, 4);
    n_array(0, 0) = ty1; n_array(0, 1) = ty2; n_array(0, 2) = ty3; n_array(0, 3) = ty4;

    return perlin1D(n_array, 0, h, (nx - i) * h + h);
}

// creating an arrary of perlin noise with multiple octaves (layers added on top of each other)
Array <float> PerlinArray(float2 size, int octaves, float scale_, float persistence, float luclarity, float min_amp, float max_amp)
{
    // the output list/scaler field
    Array <float> list = Array <float> (size.x, size.y);

    // the amplitude (not the final amplitude)
    float amplitude = 1.;

    // creating the noise (with multiple octaves)
    for (int o = 0; o < octaves; o++)
    {
        // the size for the array of random numbers
        float2 rand_size = ceil(size / scale_) + 2;
        // creating an array of random values (between 1 and 0)
        Array <float> rNoise = array::random <float> (rand_size.x, rand_size.y, 0, 1);
        // creating the noise
        for (int x = 0; x < size.x; x++)
        {
            for (int y = 0; y < size.y; y++)
            {
                // getting the perlin noise, scalling it to the amplitude, and adding it to the current height
                list(x, y) += (perlin2D(rNoise, scale_, x, y) * 2. - 1.) * amplitude;
            }
        }
        
        // changing the scale and amplitude of the noise
        scale_ *= persistence;
        amplitude *= luclarity;
    }

    // finding the max and min values of the array
    float height;
    float min_ = 999999999;
    float max_ = -999999999;
    // looping through the x and y on the image
    for (int x = 0; x < size.x; x++)
    {
        for (int y = 0; y < size.y; y++)
        {
            // finding if the height here is greater or less then the max and min
            height = list(x, y);
            min_ = (height < min_) ? height: min_;
            max_ = (height > max_) ? height: max_;
        }
    }

    // remapping the heights to your desired mins and maxes
    float new_height;
    for (int x = 0; x < size.x; x++)
    {
        for (int y = 0; y < size.y; y++)
        {
            new_height = list(x, y) - min_;
            float divisor = (max_ - min_);
            // preventing nan from popping up from division by zero
            if (divisor != 0)
            {
                new_height *= (max_amp - min_amp) / divisor;
                new_height += min_amp;
            }
            else new_height = min_amp;
            list(x, y) = new_height;
        }
    }

    // returning the heights
    return list;
}


int main()
{
    // setting up a random number generator
    srand((unsigned)time(NULL));
    
    int num_drops = 1500000. * 7.7160494;  // the number of drops to be simulated (you can both under and over erode)

    // the size of the world being simulated
    float2 world_size = float2(1500, 1500);
    // the final image (will be outputed as a custom file type to be converted to a png in python3)
    Image::Image image = Image::Image(world_size);

    // getting the perlin noise array
    Array <float> height_map = PerlinArray(world_size, 6, 250., 0.5, 0.45, 0., 1.);
    
    // so the drops wont try to access values outside of the array (by 1)
    float2 map_size = world_size - 1.;

    // the settings (i found these settings work pretty well)
    float dt = 1.2;  // simulation detail/delta time  smaller = more accurate but slower
    float density = 1.;  // the density of the water dropplet
    float minVol = 0.01;  // the minimum volume before the drop is considered to be fully evapoated
    float evapRate = 0.001;  // the rate at which the drop evaporates
    float friction = 0.05;  // the friction on the drop
    float depositionRate = 0.1;  // the rate at which the drop deposites and collects sediment
    
    // precalculating some values
    float scaled_evap_rate = (1. - dt * evapRate);
    float scaled_deposition_rate = dt * depositionRate;
    float friction_coefficient = (1. - dt * friction);

    // simulating the drops
    for (int drop_num = 0; drop_num < num_drops; drop_num++)
    {
        // finding a random position and initial speed/velocity
        float2 random_pos = float2(((float) rand() / RAND_MAX) * map_size.x, ((float) rand() / RAND_MAX) * map_size.y);
        float2 random_speed = float2(((float) rand() / RAND_MAX) * 3. - 1.5, ((float) rand() / RAND_MAX) * 3. - 1.5);
        // creating the drop
        Particle drop = Particle(random_pos, random_speed);
        
        // some values using in the simulation
        float3 n;
        float2 ipos;
        float sdiff;
        float maxsediment;

        // simulating the drop
        while (drop.volume > minVol)
        {
            ipos = floor(drop.pos);  // the grid/array position
            n = SurfaceNormal(ipos.x, ipos.y, map_size, height_map);  // finding the surface normal at this position

            drop.speed += float2(n.x, n.z) * (1. / (drop.volume * density)) * dt;  // accellerating the drop
            drop.pos   += drop.speed * dt;  // moving the drops position
            drop.speed *= friction_coefficient;  // adding friction
            
            // checking if the drop has left the map
            if (drop.pos.x < 0. || drop.pos.x > map_size.x || drop.pos.y < 0. || drop.pos.y > map_size.y) break;
            
            // not sure entirely what these next 3 lines do but i think its finding the max amount of soil the drop can collect/deposite
            maxsediment = drop.volume * length(drop.speed) * (height_map(ipos.x, ipos.y) - height_map(drop.pos.x, drop.pos.y));

            if (maxsediment < 0.) maxsediment = 0.;

            sdiff = maxsediment - drop.sediment;

            // changes the drops current sediment
            drop.sediment += scaled_deposition_rate * sdiff;
            // erroding the map / displacing the sediment on the map/height map
            height_map(ipos.x, ipos.y) -= dt * drop.volume * depositionRate * sdiff;
            
            // evaporating the drop
            drop.volume *= scaled_evap_rate;
        }
    }

    // filling out the image
    for (int x = 0; x < world_size.x; x++)
    {
        for (int y = 0; y < world_size.y; y++)
        {
            image(float2(x, y)) = float3(height_map(x, y));
        }
    }

    // saving the image as a file
    image.write("Erroded_map");  // the name of the file
}

