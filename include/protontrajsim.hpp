/*
  File contains code for proton trajectory simulation. Goal is to find correction
  of target magnet's B0 to the recoil proton angle.
 */

#pragma once

#include <algorithm>

#include <TVector3.h>
#include "srccommon.hpp"


namespace s13 {
    namespace sim {

        /// global constants
        constexpr Double_t gk_distance_step = 5e-5; // m
        constexpr Double_t gk_speed_of_light = 3e8; // m/s
        constexpr Double_t gk_momentum_mevc_to_si = 5.36e-22;
        constexpr Double_t gk_e_charge = 1.60217662e-19; // Couloumb
        constexpr Double_t gk_mass_mev_to_kg = 1.79e-30;

        class Particle {
        public:

            /*
             * Constructs new particle.
             *
             * Params:
             *   mass_mev: Mass in MeV
             *   charge: Charge in units of electron charge_e
             *   tke: Total kinetic energy of the particle in MeV
             *   name: name of the constructed particle
             */
          Particle(Double_t mass_mev, Double_t charge, std::string name)  :
                    m_name{name}, m_mass_mev{mass_mev}, m_charge{charge}, m_position{0,0,0},
                    m_direction{0,0,1} {}

            const std::string& name() const {
                return m_name;
            }

            /*
             * Returns total kinetic energy of the particle
             */
            Double_t tke_mev() const {
                return (gamma() - 1) * mass_mev();
            }

            void set_tke_mev(Double_t new_tke) {
                Double_t gamma = new_tke / mass_mev() + 1.;
                m_beta = std::sqrt(1. - (1. / std::pow(gamma, 2)));
            }

            /*
             * Returns Lorentz's gamma for the particle
             */
            Double_t gamma() const {
                return 1. / (std::sqrt(1. - std::pow(m_beta, 2)));
            }

            Double_t beta() const {
                return m_beta;
            }

            /*
             * Returns particle velocity (a scalar) in m/s
             */
            Double_t velocity() const {
                return beta() * gk_speed_of_light;
            }

            /*
             * Returns particle speed vector
             */
            TVector3 speed() const {
                return velocity() * direction();
            }

            TVector3 momentum_mevc() const {
              return gamma() * mass_mev() * beta() * direction();
            }

            TVector3 momentum_si() const {
              return momentum_mevc() * gk_momentum_mevc_to_si;
            }

            /*
             * Sets new momentum_mevc vector for the particle.
             * Updates rest of variables.
             */
            void set_momentum_mevc(TVector3 new_momentum) {
                m_direction = new_momentum * (1./new_momentum.Mag());
                m_beta = new_momentum.Mag() / std::sqrt(std::pow(mass_mev(), 2) + new_momentum.Mag2());
            }

            /*
             * Returns mass of the particle in MeV
             */
            Double_t mass_mev() const {
                return m_mass_mev;
            }

            /*
             * Returns mass of the particle in kg
             */
            Double_t mass() const {
                return m_mass_mev * gk_mass_mev_to_kg;
            }

            Double_t charge_e() const {
                return m_charge;
            }

            Double_t charge_couloumb() const {
                return m_charge * gk_e_charge;
            }

            const TVector3& position() const {
                return m_position;
            }

            void transform_position(const TVector3& dr) {
                m_position += dr;
            }

            void set_position(const TVector3& new_position) {
                m_position = new_position;
            }

            const TVector3& direction() const {
                return m_direction;
            }

            void transform_direction(const TVector3& dd) {
                m_direction += dd;
            }

            void set_direction(const TVector3& new_direction) {
                m_direction = new_direction;
            }

        private:
            std::string m_name;
            Double_t m_mass_mev;
            Double_t m_charge;
            TVector3 m_position;
            TVector3 m_direction;

            Double_t m_beta;
        };

        std::ostream& operator<<(std::ostream& ostr, const Particle& rhs);

        class Material {
        public:

            Material(std::string name, std::unique_ptr<ROOT::Math::Interpolator>&& stopping_power) :
              m_name{name}, m_stopping_power{std::move(stopping_power)}
            {}

            /*
             * Computes dE in this material given initial energy of a proton.
             *
             * Params:
             *   e_initial: Initial value of proton energy in MeV.
             *   distance: How much a particle should travel in mm.
             * Returns: dE value in MeV
             */
            virtual inline Double_t compute_dE(Double_t e_initial, Double_t distance, Double_t step_size = 0.1) {
                Double_t cur_dist = 0;
                Double_t cur_e = e_initial;
                Double_t de = 0;
                while(cur_dist < distance) {
                  if (cur_e < 1.1) {
                    de = e_initial;
                    return de;
                  }

                  Double_t step_de = m_stopping_power->Eval(cur_e) * step_size;
                  cur_e -= step_de;
                  de += step_de;

                  cur_dist += step_size;
                }

                if (de > e_initial) {
                  de = e_initial;
                }
                return de;
            }

            const std::string& name() const {
                return m_name;
            }

            bool operator==(const std::string& other_name);

        protected:
            std::string m_name;
            std::shared_ptr<ROOT::Math::Interpolator> m_stopping_power;
            Double_t m_step_size;
        };

        bool operator==(const Material& lhs, const Material& rhs);

        class MaterialManager {
        public:
            MaterialManager() {}

            void Add(std::shared_ptr<Material> material) {
                m_materials.push_back(material);
            }

            std::shared_ptr<Material> Find(std::string name) {
                auto it = std::find_if(m_materials.begin(), m_materials.end(),
                                       [&name] (std::shared_ptr<Material> mtrl) -> bool {
                    return mtrl->name() == name;
                });
                if (it == m_materials.end()) {
                    throw std::invalid_argument("Material with the given name does not exists!");
                }

                return *it;
            }
        private:
            std::vector<std::shared_ptr<Material> > m_materials;
        };

        std::shared_ptr<Material> LoadMaterial(std::string name, std::string stp_power_file);

        extern s13::sim::MaterialManager g_material_manager;
        extern std::shared_ptr<s13::sim::Material> g_default_material;

        enum class ZOrderPolicy {
            Exclude,
            Overlay
        };

      class Volume {
      public:

        Volume(std::string name) : m_name{name}, m_material{g_default_material},
                                   m_zorder_policy{ZOrderPolicy::Exclude}{}

        const std::string& name() {
          return m_name;
        }

        /*
          Tells if point is inside the object's boundaries.
        */
        virtual bool Contains(const TVector3& pos) const = 0;

        /*
         * Computes dE for the current simulation step.
         *
         * Params:
         *   pos: Vector of proton position to test if it is in this geometry
         *   e_initial: Initial value of proton energy
         *   step_size: size of the step on which dE will be computed
         * Returns:
         *   Value of dE on this step if the proton is inside this geometry otherwise
         *   returns 0.
         */
        Double_t ComputeDe(const TVector3& pos, Double_t e_initial) {
          Double_t dist_step_mm = gk_distance_step * 1e3;
          if (Contains(pos)) {
            double de = m_material->compute_dE(e_initial, dist_step_mm,
                                               dist_step_mm / 3);
            // if (m_name == "chamber") {
            //   std::cout << "dE in chamber: " << de << std::endl;
            //   std::cout << "Input args to Material::compute_dE(" << e_initial
            //             << ", " << dist_step_mm << ", " << dist_step_mm / 3
            //             << ");" << std::endl;
            //   std::cout << "Material name: " << m_material->name() << std::endl;
            //   std:: cout << "Stopping power at this e_initial: "
            //              << m_material->m_stopping_power->Eval(e_initial) << std::endl;
            // }
            m_tot_de += de;
            return de;
          } else {
            if (m_name == "chamber") {
              std::cout << "returning 0 dE in chamber" << std::endl;
            }
            return 0;
          }
        }

          /*
           * Updates given momentum_mevc vector of the proton.
           *
           * This function should account for all effects of the object, like dE,
           * magnetic field (if present), etc.
           */
          virtual void UpdateDirection(Particle& particle) {
            // empty implementation for default case
          }

          virtual void UpdateEnergy(Particle& particle) {
            Double_t de = ComputeDe(particle.position(), particle.tke_mev());
            particle.set_tke_mev(particle.tke_mev() - de);
          }

          bool Transform(Particle& particle) {
            if (!Contains(particle.position())) {
              // transform was not applied
              return false;
            }
            UpdateDirection(particle);
            UpdateEnergy(particle);
            return true;
          }

          ZOrderPolicy GetZOrderPolicy() {
            return m_zorder_policy;
          }

          std::shared_ptr<Material> GetMaterial() {
            return m_material;
          }

        void SetMaterial(std::shared_ptr<Material> material) {
          m_material = material;
        }

        double GetTotalDe() {
          return m_tot_de;
        }

        /*
          Resets all accumulating metrics of a volume
        */
        void Reset() {
          m_tot_de = 0;
        }

        protected:
          std::string m_name;
          std::shared_ptr<Material> m_material;
          double m_tot_de = 0;
          ZOrderPolicy m_zorder_policy;
        };

        class World : public Volume {
        public:
            World(std::shared_ptr<Material> material = g_default_material) : Volume("WORLD") {
                m_material = material;
            };

            virtual bool Contains(const TVector3 &pos) const {
                return true;
            }
        };

        class TargetCrystal : public Volume {
        public:
            TargetCrystal(std::string name, Double_t radius, Double_t thickness, TVector3 center_pos,
                          std::shared_ptr<Material> material = g_default_material) :
                    Volume{name}, m_radius{radius}, m_thickness{thickness}, m_position{center_pos} {
                m_material = material;
            }

            virtual bool Contains(const TVector3 &pos) const {
                Double_t obj_r = std::sqrt(std::pow(pos.X() - m_position.X(), 2) +
                                           std::pow(pos.Y() - m_position.Y(), 2));
                return obj_r <= m_radius &&
                       std::abs(pos.Z() - m_position.Z()) < m_thickness / 2;
            }

        private:
            Double_t m_radius;
            Double_t m_thickness;
            TVector3 m_position;
        };


        /*
         * Defines horizontally oriented cylinder volume
         */
        class HorizontalCylinder : public Volume {
        public:
            HorizontalCylinder(std::string name, Double_t radius, Double_t height, TVector3 center_pos,
                               std::shared_ptr<Material> material = g_default_material) :
              Volume{name}, m_radius{radius}, m_height{height}, m_position{center_pos} {
                m_material = material;
                m_zorder_policy = ZOrderPolicy::Exclude;
            };

            virtual bool Contains(const TVector3 &pos) const {
                Double_t obj_r = std::sqrt(std::pow(pos.X() - m_position.X(), 2) +
                                           std::pow(pos.Z() - m_position.Z(), 2));
                bool contains = obj_r <= m_radius &&
                  std::abs(pos.Y() - m_position.Y()) < m_height / 2;

                return contains;
            }

        private:
            Double_t m_radius;
            Double_t m_height;
            TVector3 m_position;
        };

        /*
         * Defines magnetic field in a cylindrical volume.
         */
        class HorizontalCylinderWithMagField : public HorizontalCylinder {
        public:
            HorizontalCylinderWithMagField(std::string name, Double_t radius, Double_t height,
                                           TVector3 center_pos, TVector3 b0,
                                           std::shared_ptr<Material> material = g_default_material) :
              HorizontalCylinder{name, radius, height, center_pos, material}, m_b0{b0} {
                m_zorder_policy = ZOrderPolicy::Overlay;
            };

            virtual void UpdateDirection(Particle& particle) {
              if (particle.velocity() < 0.1) {
                return;
              }
                Double_t dt = gk_distance_step / particle.velocity();
                // std::cout << "dt: " << dt << " s" << std::endl;

                TVector3 b0_force = particle.charge_couloumb() * particle.speed().Cross(m_b0);
//                std::cout << "particle speed: " << particle.speed().Mag() << " m/s" << std::endl;
//                std::cout << "particle charge: " << particle.charge_couloumb() << " C" << std::endl;

                TVector3 b0_force_momentum_mevc = (b0_force * dt) * (1 / gk_momentum_mevc_to_si);
//                std::cout << "b0_force_momentum_mevc: " << b0_force_momentum_mevc.Mag() << " MeV/c" << std::endl;

                TVector3 new_particle_momentum_mevc = particle.momentum_mevc() + b0_force_momentum_mevc;
//                std::cout << "particle theta: " <<
//                          new_particle_momentum_mevc.Theta() * 1/s13::common::D2R << std::endl;
                particle.set_momentum_mevc(new_particle_momentum_mevc);
            }

            TVector3 B0() {
                return m_b0;
            }

            void B0(TVector3 new_val) {
                m_b0 = new_val;
            }

        private:
            Double_t m_radius;
            Double_t m_height;
            TVector3 m_position;
            TVector3 m_b0;
        };

        /*
         * Defines RDC box to detect proton hitting it.
         */
        class RectangularBox : public Volume {
        public:
            RectangularBox(std::string name,
                   Double_t widthX, Double_t heightY, Double_t depthZ,
                   Double_t dist_from_origin, Double_t rotationY) :
                    Volume(name),
                    m_widthX{widthX}, m_heightY{heightY}, m_depthZ{depthZ},
                    m_dist_from_origin{dist_from_origin}, m_rotationY{rotationY} {}

            virtual bool Contains(const TVector3 &pos) const {
                TVector3 point = pos;
                // transform given point into local coordinate frame
                point.RotateY(-m_rotationY);
                point.SetZ(point.Z() - m_dist_from_origin);

                if (std::abs(point.X()) < m_widthX / 2 &&
                        std::abs(point.Y()) < m_heightY /2 &&
                        std::abs(point.Z()) < m_depthZ / 2) {
                    return true;
                }

                return false;
            }

          TVector3 to_lcs(TVector3 point) const {
            point.RotateY(-m_rotationY);
            point.SetZ(point.Z() - m_dist_from_origin);
            return point;
          }

        private:
            Double_t m_widthX;
            Double_t m_heightY;
            Double_t m_depthZ;
            Double_t m_dist_from_origin;
            Double_t m_rotationY;
        };

        class RunManager {
        public:

            RunManager() {}

            void AddGeometry(std::shared_ptr<Volume> geometry) {
                m_geometries.push_back(geometry);
            }

            void AddParticle(std::shared_ptr<Particle> particle) {
                m_particles.push_back(particle);
            }

            void ClearParticles() {
                m_particles.clear();
            }

            std::shared_ptr<Volume> FindGeometry(std::string name) {
                auto it = std::find_if(m_geometries.begin(), m_geometries.end(),
                                       [&name] (std::shared_ptr<Volume> geom) -> bool {
                                           return geom->name() == name;
                                       });
                if (it == m_geometries.end()) {
                    throw std::invalid_argument("Volume with the given name does not exists!");
                }

                return *it;
            }

            std::shared_ptr<Particle> FindParticle(std::string name) {
                auto it = std::find_if(m_particles.begin(), m_particles.end(),
                                       [&name] (std::shared_ptr<Particle> part) -> bool {
                                           return part->name() == name;
                                       });
                if (it == m_particles.end()) {
                    throw std::invalid_argument("Particle with the given name does not exists!");
                }

                return *it;
            }

            template <typename UnaryFn>
            void Iterate(UnaryFn fn) {
                for (auto iterator = m_geometries.begin(); iterator != m_geometries.end(); ++iterator) {
                    bool will_continue = fn(*iterator);
                    if (!will_continue) {
                        break;
                    }
                }
            }

            void ProcessStep() {
                ProcessGeometries();
                UpdateParticlePositions();
            }

          void ResetVolumes() {
            for (auto v : m_geometries) {
              v->Reset();
            }
          }

        protected:

            void UpdateParticlePositions() {
                for (auto iterator = m_particles.begin(); iterator != m_particles.end(); ++iterator) {
                    auto particle= *iterator;
                    TVector3 dr = gk_distance_step * particle->direction();
                    particle->transform_position(dr);
                }
            }

            void ProcessGeometries() {
                bool transform_applied = false;
                Iterate([&] (std::shared_ptr<Volume> volume) -> bool {
                    for (auto iterator = m_particles.begin(); iterator != m_particles.end(); ++iterator) {
                        auto particle = *iterator;
                        if (volume->GetZOrderPolicy() == ZOrderPolicy::Exclude && !transform_applied) {
                            // if it is solid, gas - a substance, then only a Volume first in the z-order
                            // should apply its transforms to a particle
                            transform_applied = volume->Transform(*particle);
                        } else if (volume->GetZOrderPolicy() == ZOrderPolicy::Overlay) {
                            // if it is a field, i.e. it "overlays" other Volumes, then apply its transform
                            // regardless of how many previous transforms were already applied
                          volume->Transform(*particle);
                        }
                    }
                    return true;
                });
            }

        private:
            std::vector<std::shared_ptr<Volume> > m_geometries;
            std::vector<std::shared_ptr<Particle> > m_particles;
        };

        extern s13::sim::RunManager g_run_manager;

    /*
      Computes angle as measured by RDC by given angle (from actual kinematics).
      Function goal is to account for proton trajectory distortion by the
      target magnet's B0 and e-loss in different materials.

      Params:
        proton_theta: Recoil proton scattering angle in radians.
        vertex_pos: Vector with initial coordinates of the proton.
      Returns:
        Proton angle as it is measured by the RDC.
    */
      std::vector<Double_t> proton_theta_measured(std::vector<Double_t> proton_theta, TVector3 vertex_pos);

      std::vector<Double_t>
      proton_delta_x_vs_phi(std::vector<Double_t> proton_phi_range,
                            Double_t proton_theta,
                            TVector3 vertex_pos);
    }
}
