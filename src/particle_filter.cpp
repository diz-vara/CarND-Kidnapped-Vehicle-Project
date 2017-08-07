/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 *
 *  Solution: Aug 6, 2017 - Anton Varfolomeev
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;


Particle::Particle(int _id, double _x, double _y, double _theta, double _weight) :
  id(_id),
  x(_x),
  y(_y),
  theta(_theta), 
  weight(_weight)
{

}


void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
  if (num_particles < 1)
    std::cerr << "Invalid number of particles (" << num_particles << ")" << std::endl;


  //create independent distributions for x, y and theta
  std::normal_distribution<> d_x(x, std[0]);
  std::normal_distribution<> d_y(y, std[1]);
  std::normal_distribution<> d_t(theta, std[2]);
  for (int i = 0; i < num_particles; ++i) {
    particles.push_back(Particle(i+1,d_x(gen), d_y(gen), d_t(gen)));
    weights.push_back(1.);
  }
  is_initialized = true;


}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

  //I create three zero-centerd distributions - and will shift them to individual means
  std::normal_distribution<> d_x(0, std_pos[0]);
  std::normal_distribution<> d_y(0, std_pos[1]);
  std::normal_distribution<> d_t(0, std_pos[2]);

  for (Particle &p: particles) {
    if (yaw_rate != 0) {
      double new_theta = p.theta + yaw_rate * delta_t;
      p.x = p.x + velocity / yaw_rate * (sin(new_theta) - sin(p.theta)) + d_x(gen);
      p.y = p.y + velocity / yaw_rate * (cos(p.theta) - cos(new_theta)) + d_y(gen);
      p.theta = new_theta + d_t(gen);
    }
    else {
      p.x = p.x + velocity * cos(p.theta) * delta_t + d_x(gen);
      p.y = p.y + velocity * sin(p.theta) * delta_t + d_y(gen);
      p.theta = p.theta + d_t(gen);
    }
  }

}

// 2-d normal distribution probability
double ParticleFilter::Mv_NormP(double x, double y, double xm, double ym, double s_x, double s_y) {
  double s1 = (x - xm)*(x - xm) / (2. * s_x * s_x);
  double s2 = (y - ym)*(y - ym) / (2. * s_y * s_y);
  double e = exp(-1. * (s1 + s2));
  double w = e / (2. * M_PI * s_x * s_y);
  return w;
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	double s_x = std_landmark[0];   //sigma x
	double s_y = std_landmark[1];	//sigma y
	std::normal_distribution<> d_x(0, s_x);
	std::normal_distribution<> d_y(0, s_y);

  int i(0);
  for ( Particle &p: particles) {
    
    std::vector<int> associations;
    std::vector<double> sense_x;
    std::vector<double> sense_y;
    double weight = 1;
    double cos_t(cos(p.theta));
		double sin_t(sin(p.theta));

    for (LandmarkObs obs : observations) {
			//translate from vehicle to Map coordinates
			double x = obs.x * cos_t - obs.y * sin_t + p.x + d_x(gen); 
			double y = obs.x * sin_t + obs.y * cos_t + p.y + d_y(gen);
			double min_dist(sensor_range), cur_dist;
			int idx(-1);
			Map::single_landmark_s best_landmark;
			for (Map::single_landmark_s lmk : map_landmarks.landmark_list) {
        cur_dist = dist(x, y, lmk.x_f, lmk.y_f);
				if (cur_dist < min_dist) {
					min_dist = cur_dist;
					idx = lmk.id_i;
					best_landmark = lmk;
				}
			}
			if (idx >= 0 ) {
				associations.push_back(idx);
				sense_x.push_back(x);
				sense_y.push_back(y);
        weight *= Mv_NormP(x, y, best_landmark.x_f, best_landmark.y_f, s_x, s_y);
			}
		}

    //no valid associations!!!
    if (associations.empty())
      weight = 0;

    weights[i++] = weight;
    p = SetAssociations(p, associations, sense_x, sense_y);
    p.weight = weight;
  }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	std::discrete_distribution<unsigned int> d(weights.begin(), weights.end());

	std::vector<Particle> oldParticles(particles);
	for (int i = 0; i < num_particles; ++i) {
		particles[i] = oldParticles[d(gen)];
  }
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
