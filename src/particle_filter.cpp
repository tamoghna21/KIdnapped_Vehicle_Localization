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

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles =10;
	default_random_engine gen;
	
	// creates a normal (Gaussian) distribution for x, y, theta
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);
	
	for (unsigned int i =0; i< num_particles; i++){
	    Particle particle;
	    
		particle.x = dist_x(gen);
		particle.y = dist_y(gen);
		particle.theta = dist_theta(gen);
		
		particles.push_back(particle);
    }
	
    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	default_random_engine gen;
	Particle particle;
	double xf;
	double yf;
	double thetaf;
	
	for (unsigned int i =0; i< num_particles; i++){
	    particle = particles[i];
	    
	    if(fabs(yaw_rate) < 0.0001){
	    	xf = particle.x + velocity * delta_t * cos(particle.theta);
			yf = particle.y + velocity * delta_t * sin(particle.theta);
			thetaf = particle.theta;
		}
		else{
			xf = particle.x + (velocity/yaw_rate)*(sin(particle.theta+(yaw_rate*delta_t))-sin(particle.theta));
			yf = particle.y + (velocity/yaw_rate)*(cos(particle.theta) - cos(particle.theta+(yaw_rate*delta_t)));
			thetaf = particle.theta + yaw_rate * delta_t;
		}
		
		normal_distribution<double> dist_x(xf, std_pos[0]);
		normal_distribution<double> dist_y(yf, std_pos[1]);
		normal_distribution<double> dist_theta(thetaf, std_pos[2]);
		
		particles[i].x = dist_x(gen);
		particles[i].y = dist_y(gen);
		particles[i].theta = dist_theta(gen); 
		
		//cout << "particle x: " << particles[i].x << endl;
		//cout << "particle y: " << particles[i].y << endl;
		//cout << "particle theta: " << particles[i].theta << endl;
		
    }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> &predicted, const std::vector<LandmarkObs>& trans_observations,const Map &map_landmarks) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	double dist;
	double min_dist;
	LandmarkObs predicted_landmark;
	for (unsigned int i=0; i< trans_observations.size(); i++){
	    min_dist = 1000000.0;
	    predicted_landmark.x = 0.;
	    predicted_landmark.y = 0.;
		for (unsigned int k=0; k< map_landmarks.landmark_list.size(); k++){
			dist = sqrt(pow((trans_observations[i].x - map_landmarks.landmark_list[k].x_f),2.0)+pow((trans_observations[i].y - map_landmarks.landmark_list[k].y_f),2.0));
		
			if (dist < min_dist){
				min_dist = dist;
				predicted_landmark.id = map_landmarks.landmark_list[k].id_i;
				predicted_landmark.x = map_landmarks.landmark_list[k].x_f;
				predicted_landmark.y = map_landmarks.landmark_list[k].y_f;
			}
		}
		predicted.push_back(predicted_landmark);
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
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
	
	std::vector<LandmarkObs> predicted_landmarks;
	std::vector<LandmarkObs> trans_observations;
	
	LandmarkObs s_trans_obs;
	double mul_prob;
	double top1;
	double top2;
	//double xm;
	//double ym;
	
	for (unsigned int i=0; i< num_particles; i++){
		trans_observations.clear();
		predicted_landmarks.clear();
		mul_prob = 1.0;
		for (unsigned int j=0; j< observations.size(); j++){
			// Homogeneous transformation
			s_trans_obs.x  = ((observations[j].x) * cos(particles[i].theta)) - ((observations[j].y) * sin(particles[i].theta)) + particles[i].x;
			s_trans_obs.y  = ((observations[j].x) * sin(particles[i].theta)) + ((observations[j].y) * cos(particles[i].theta)) + particles[i].y;
			
			trans_observations.push_back(s_trans_obs);
		}
		//Data association
		dataAssociation(predicted_landmarks, trans_observations,map_landmarks);
		
		//particles[i] = SetAssociations(particles[i], const std::vector<int>& associations, 
        //                             const std::vector<double>& sense_x, const std::vector<double>& sense_y)
		
		for (unsigned int j=0; j< trans_observations.size(); j++){
		    top1 = (pow(trans_observations[j].x - predicted_landmarks[j].x, 2.0))/ (2*std_landmark[0]*std_landmark[0]);
		    top2 = (pow(trans_observations[j].y - predicted_landmarks[j].y, 2.0))/ (2*std_landmark[1]*std_landmark[1]);
			mul_prob *= (1 / (2 * M_PI * std_landmark[0] * std_landmark[1]))* exp(-1.0 * (top1+top2));
		}
		
		particles[i].weight = mul_prob;
	}
	
	
		
	
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	
	std::vector<double> weights;
	std::vector<Particle> particles_sampled;
	default_random_engine gen;
	//std::discrete_distribution<int> distribution {2,2,1,1,2,2,1,1,2,2};
	
	for (unsigned int i=0; i< num_particles; i++){
		weights.push_back(particles[i].weight);
	}
	discrete_distribution<> distribution(weights.begin(), weights.end());
	
	for (unsigned int i=0; i< num_particles; ++i){
		int number = distribution(gen);
		particles_sampled.push_back(particles[number]);
	}
	particles = particles_sampled;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
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
