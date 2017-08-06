/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
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

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	for (LandmarkObs obs : observations) {
		double min_dist = 1e9; //large enough value
		for (LandmarkObs prd : predicted) {

		}
	}

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
    p.associations.clear();
    p.sense_x.clear();
    p.sense_y.clear();
    double weight = 1;
    double cos_t(cos(p.theta));
		double sin_t(sin(p.theta));
		// go through all observations
    bool bObserved(false);

    for (LandmarkObs obs : observations) {
			//transfor from vehicle to Map coordinates
			double x = obs.x * cos_t - obs.y * sin_t + p.x + d_x(gen); //todo: noise
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
        bObserved = true;
				p.associations.push_back(idx);
				p.sense_x.push_back(x);
				p.sense_y.push_back(y);
        double s1 = (x - best_landmark.x_f)*(x - best_landmark.x_f) / (2. * s_x * s_x);
        double s2 = (y - best_landmark.y_f)*(y - best_landmark.y_f) / (2. * s_y * s_y);
        double pow = -1. * (s1 + s2);
        double e = exp(pow);
        double w = e / (2. * M_PI * s_x * s_y);
        weight *= w;
			}
		}

    if (!bObserved)
      weight = 0;
    weights[i++] = weight;
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
