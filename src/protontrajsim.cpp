
#include <vector>
#include <memory>

#include "protontrajsim.hpp"


using namespace s13::common;
using namespace s13::sim;


s13::sim::MaterialManager s13::sim::g_material_manager;
std::shared_ptr<s13::sim::Material> s13::sim::g_default_material = LoadMaterial("vacuum", "data/p_in_vacuum.csv");
s13::sim::RunManager s13::sim::g_run_manager;


bool s13::sim::operator==(const Material &lhs, const Material &rhs) {
    return lhs.name() == rhs.name();
}

bool s13::sim::Material::operator==(const std::string &other_name) {
    return other_name == m_name;
}

std::shared_ptr<Material> s13::sim::LoadMaterial(std::string name, std::string stp_power_file) {
    auto stp_power = s13::common::LoadStoppingPower(stp_power_file);
    std::shared_ptr<Material> material(new Material(name, std::move(stp_power)));
    g_material_manager.Add(material);
    return material;
}


///////////////////////
/// Simulation code
///////////////////////

enum class Espri {
    left,
    right
};

enum class ReactionType {
    p_6he,
    p_4he
};

static Espri g_which_espri = Espri::left;
static ReactionType g_reaction_type = ReactionType::p_6he;

// for theta_eff determination & losses in target materials
// void SetupGeometry() {
//     Double_t espri_center_angle = 62.5;
//     Double_t dist_to_rdc_plane = 1005e-3;
//     if (g_which_espri == Espri::right) {
//         espri_center_angle = -espri_center_angle;
//     }

//     std::shared_ptr<Volume> world{new World(LoadMaterial("air", "data/p_in_air.csv"))};
//     std::shared_ptr<Volume> target{
//       new TargetCrystal("target", 12e-3, 3e-3, TVector3(0, 0, 0),
//                         LoadMaterial("naphthahlene", "data/p_in_naph.csv"))};
//     std::shared_ptr<Volume> chamber{
//       new HorizontalCylinder("chamber", 109e-3, 80e-3, TVector3(0, 0, 0),
//                              LoadMaterial("N2", "data/p_in_N2_gas.csv"))};
//     std::shared_ptr<Volume> magnet{
//       new HorizontalCylinderWithMagField("magnet", 110e-3, 90e-3,
//                                          TVector3(0, 0, 0),
//                                          TVector3(0,-96.12e-3,0))};
//     std::shared_ptr<Volume> nai_crystals_plane{
//       new RectangularBox("RDCBox", 840e-3, 840e-3, 1e-3, dist_to_rdc_plane,
//                  espri_center_angle * D2R)};

//     g_run_manager.AddGeometry(target);
//     g_run_manager.AddGeometry(chamber);
//     g_run_manager.AddGeometry(magnet);
//     g_run_manager.AddGeometry(nai_crystals_plane);
//     g_run_manager.AddGeometry(world);
// }

// to compute E spectra on NaIs including all e-losses
// additionally to SetupGeometry it includes 4 mm plastic
// scintillator and 1 mm Al cover of NaIs
void SetupGeometry2() {
    Double_t espri_center_angle = 62.5;
    Double_t dist_to_rdc_plane = 1005e-3;
    Double_t dist_to_espri_pla = dist_to_rdc_plane + 100e-3;
    Double_t dist_to_nais_plane = dist_to_rdc_plane + 340e-3;
    Double_t dist_to_nais_crystals = dist_to_nais_plane + 1e-3;
    if (g_which_espri == Espri::right) {
        espri_center_angle = -espri_center_angle;
    }

    std::shared_ptr<Volume> world{new World(LoadMaterial("air", "data/p_in_air.csv"))};
    std::shared_ptr<Volume> target{
      new TargetCrystal("target", 12e-3, 3e-3, TVector3(0, 0, 0),
                        LoadMaterial("naphthahlene", "data/p_in_naph.csv"))};
    std::shared_ptr<Volume> chamber{
      new HorizontalCylinder("chamber", 109e-3, 80e-3, TVector3(0, 0, 0),
                             LoadMaterial("N2", "data/p_in_N2_gas.csv"))};
    std::shared_ptr<Volume> magnet{
      new HorizontalCylinderWithMagField("magnet", 110e-3, 90e-3,
                                         TVector3(0, 0, 0),
                                         TVector3(0,-96.12e-3,0))};
    std::shared_ptr<Volume> espri_plastic_box{
      new RectangularBox("ESPRIPla", 840e-3, 840e-3, 4e-3, dist_to_espri_pla,
                         espri_center_angle * D2R)};
    espri_plastic_box->SetMaterial(LoadMaterial("BC400", "data/p_in_BC400.csv"));
    std::shared_ptr<Volume> nai_box{
      new RectangularBox("NaiAlCover", 840e-3, 840e-3, 0.11e-3, dist_to_nais_plane,
                         espri_center_angle * D2R)};
    nai_box->SetMaterial(LoadMaterial("Al", "data/p_in_Al.csv"));
    std::shared_ptr<Volume> nai_crystals_plane{
      new RectangularBox("NaiCrystalsPlane", 840e-3, 840e-3, 1e-3, dist_to_nais_crystals,
                         espri_center_angle * D2R)};

    g_run_manager.AddGeometry(target);
    g_run_manager.AddGeometry(chamber);
    g_run_manager.AddGeometry(magnet);
    g_run_manager.AddGeometry(espri_plastic_box);
    g_run_manager.AddGeometry(nai_box);
    g_run_manager.AddGeometry(nai_crystals_plane);
    g_run_manager.AddGeometry(world);
}

void SetupGeometryProtonDeltaXGeometry() {
    Double_t espri_center_angle = 62.5;
    Double_t dist_to_rdc_plane = 1005e-3;
    Double_t dist_to_espri_pla = dist_to_rdc_plane + 100e-3;
    Double_t dist_to_nais_plane = dist_to_rdc_plane + 340e-3;
    Double_t dist_to_nais_crystals = dist_to_nais_plane + 1e-3;
    if (g_which_espri == Espri::right) {
        espri_center_angle = -espri_center_angle;
    }

    std::shared_ptr<Volume> world{new World(LoadMaterial("air", "data/p_in_air.csv"))};
    std::shared_ptr<Volume> espri_rdc{
      new RectangularBox("ESPRIRDC", 440e-3, 190e-3 * 2, 1e-3, dist_to_rdc_plane,
                         espri_center_angle * D2R)};
    espri_rdc->SetMaterial(LoadMaterial("BC400", "data/p_in_BC400.csv"));

    g_run_manager.AddGeometry(espri_rdc);
    g_run_manager.AddGeometry(world);
}

void SetupParticles(TVector3 vertex_pos) {
    Particle *proton = new Particle(938, 1, "proton");
    proton->set_position(vertex_pos);
    TVector3 new_dir = proton->direction();
    if (g_which_espri == Espri::left) {
      new_dir.RotateY(62.5 * D2R);
    } else {
        new_dir.RotateY(-62.5 * D2R);
    }
    proton->set_direction(new_dir);
    g_run_manager.AddParticle(std::shared_ptr<Particle>(proton));
}

void SetupProton(Double_t eff_angle_deg, TVector3 vertex_pos, std::shared_ptr<Particle> proton) {
    proton->set_position(vertex_pos);
    TVector3 new_dir{0,0,1};
    if (g_which_espri == Espri::left) {
        new_dir.RotateY(eff_angle_deg * D2R);
    } else {
        new_dir.RotateY(-eff_angle_deg * D2R);
    }
    proton->set_direction(new_dir);
    if (g_reaction_type == ReactionType::p_4he) {
        auto p_he4_e = LoadProtonKinE("data/he4_kin_E.csv");
        proton->set_tke_mev(p_he4_e->Eval(eff_angle_deg * D2R));
    } else {
        auto p_he6_e = LoadProtonKinE("data/he6_kin_E.csv");
        proton->set_tke_mev(p_he6_e->Eval(eff_angle_deg * D2R));
    }
    std::cout << "setting proton energy for angle: " << eff_angle_deg << "; tke: "
              << proton->tke_mev() << " MeV" << std::endl;
}

void SetupProtonPhi(Double_t eff_theta, Double_t eff_phi, TVector3 vertex_pos, std::shared_ptr<Particle> proton) {
  SetupProton(eff_theta, vertex_pos, proton);
  TVector3 dir = proton->direction();
  dir.RotateZ(eff_phi * D2R);
  proton->set_direction(dir);
}

std::ofstream SetupOutputFile() {
    std::ofstream csv_data_ofs;

    if (g_reaction_type == ReactionType::p_4he) {
        if (g_which_espri == Espri::left) {
            csv_data_ofs = std::ofstream("out/p_4he_esl_proton_traj_sim.csv");
        } else {
            csv_data_ofs = std::ofstream("out/p_4he_esr_proton_traj_sim.csv");
        }
    } else {
        if (g_which_espri == Espri::left) {
            csv_data_ofs = std::ofstream("out/p_6he_esl_proton_traj_sim.csv");
        } else {
            csv_data_ofs = std::ofstream("out/p_6he_esr_proton_traj_sim.csv");
        }
    }

    return csv_data_ofs;
}

std::vector<Double_t> s13::sim::proton_theta_measured(std::vector<Double_t> proton_phi_range,
                                                      TVector3 vertex_pos) {
    std::vector<Double_t> measured_angles_range;
    std::ofstream csv_data_ofs = SetupOutputFile();
    SetupGeometry2();
    SetupParticles(vertex_pos);
    auto magnet = g_run_manager.FindGeometry("magnet");
    auto proton = g_run_manager.FindParticle("proton");
    auto target = g_run_manager.FindGeometry("target");
    auto chamber = g_run_manager.FindGeometry("chamber");
    auto espri_pla = g_run_manager.FindGeometry("ESPRIPla");
    auto nais_al_cover = g_run_manager.FindGeometry("NaiAlCover");
    auto nai_crystals_plane = g_run_manager.FindGeometry("NaiCrystalsPlane");
    auto world = g_run_manager.FindGeometry("WORLD");

    csv_data_ofs << "tke," << "eff_theta," << "meas_theta,"
                 << "de_target," << "de_chamber," << "de_air,"
                 << "de_espri_pla," << "de_nais_al_cover,"
                 << "tke - dEs" << std::endl;
    for (auto angle_it = proton_phi_range.begin(); angle_it != proton_phi_range.end(); ++angle_it) {
        Double_t proton_eff_theta_deg = *angle_it;
        SetupProton(proton_eff_theta_deg, vertex_pos, proton);
        g_run_manager.ResetVolumes();

        std::cout << "Running simulation for proton angle: " << proton_eff_theta_deg << "...\n\n";
        int iteration = 0;
        Double_t tke_initial = proton->tke_mev();
        while (!nai_crystals_plane->Contains(proton->position()) || proton->position().Mag() > 4000e-3) {
            // if (iteration % 500 == 0 && magnet->Contains(proton->position())) {
            //   std::cout << *proton;
            // }
            if (proton->direction().Theta() * 1./D2R > 90) {
                std::cout << "Magnetic field is too strong, proton cannot escape the yoke region." << std::endl;
                break;
            }
            if (proton->tke_mev() < 0.1) {
              std::cout << "Proton did not reach the RDC!" << std::endl;
              break;
            }

            g_run_manager.ProcessStep();
            ++iteration;
        }
        measured_angles_range.push_back(proton->position().Theta() * 1/D2R);
        std::cout << *proton;

        csv_data_ofs << tke_initial << "," << proton_eff_theta_deg << ","
                     << proton->position().Theta() * 1/D2R << "," << target->GetTotalDe()
                     << "," << chamber->GetTotalDe() << "," << world->GetTotalDe()
                     << "," << espri_pla->GetTotalDe() << "," << nais_al_cover->GetTotalDe()
                     << "," << proton->tke_mev() << std::endl;
    }

    return measured_angles_range;
}

Double_t compute_edeg_thickness(Double_t proton_theta) {
  static const Double_t edeg_dist_from_target = 1123;
  static const Double_t edeg_base_angle = 53 * D2R;
  // static const Double_t edeg_finish_angle = 65 * D2R;
  static const Double_t edeg_thickness = 29;
  static const Double_t edeg_length = 255;
  static const Double_t edeg_min_thickness = 2;
  static const Double_t edeg_angle_tan = edeg_thickness / edeg_length;

  if (proton_theta < edeg_base_angle) {
    return 0;
  }

  Double_t dist_from_base = 2 * edeg_dist_from_target *
    TMath::Tan((proton_theta - edeg_base_angle) / 2);
  Double_t length2 = edeg_length - dist_from_base;
  Double_t thickness = length2 * edeg_angle_tan;
  if (thickness < edeg_min_thickness) {
    thickness = 0;
  }
  // std::cout << "Proton angle: " << proton_theta / D2R <<
  //   "; edeg thickness: " << thickness << std::endl;

  return thickness;
}

std::vector<Double_t> s13::sim::proton_delta_x_vs_phi(std::vector<Double_t> proton_phi_range,
                                                      Double_t proton_theta,
                                                      TVector3 vertex_pos) {
  std::vector<Double_t> rdc_xes;
  std::ofstream csv_data_ofs = SetupOutputFile();
  SetupGeometryProtonDeltaXGeometry();
  SetupParticles(vertex_pos);
  auto proton = g_run_manager.FindParticle("proton");
  auto espri_rdc = std::dynamic_pointer_cast<RectangularBox>
    (g_run_manager.FindGeometry("ESPRIRDC"));
  auto world = g_run_manager.FindGeometry("WORLD");

  csv_data_ofs << "tke," << "theta," << "phi," << "Xrdc" << ",Yrdc,Zedeg" << std::endl;
  for (auto angle_it = proton_phi_range.begin(); angle_it != proton_phi_range.end(); ++angle_it) {
    Double_t proton_phi = *angle_it;
    SetupProtonPhi(proton_theta, proton_phi, vertex_pos, proton);
    g_run_manager.ResetVolumes();

    std::cout << "Running simulation for proton phi angle: " << proton_phi << "...\n\n";
    std::cout << "Proton phi angle: " << proton->direction().Phi() / D2R << std::endl;
    int iteration = 0;
    Double_t tke_initial = proton->tke_mev();
    while (!espri_rdc->Contains(proton->position())) {
      // if (iteration % 500 == 0 && magnet->Contains(proton->position())) {
      //   std::cout << *proton;
      // }
      if (proton->position().Mag() > 2000e-3) {
        break;
      }
      if (proton->tke_mev() < 0.1) {
        std::cout << "Proton did not reach the RDC!" << std::endl;
        break;
      }

      g_run_manager.ProcessStep();
      ++iteration;
    }
    if (!espri_rdc->Contains(proton->position())) {
      continue;
    }
    auto proton_rdc_pos = espri_rdc->to_lcs(proton->position());
    auto aang_vec = TVector3(proton->position().X(), 0, proton->position().Z());
    Double_t edeg_thickness = compute_edeg_thickness(aang_vec.Theta());
    rdc_xes.push_back(proton_rdc_pos.X());
    std::cout << *proton;

    std::cout << "Proton X pos at RDC: " << proton_rdc_pos.X() * 1e3 << std::endl;
    std::cout << "Proton Y pos at RDC: " << proton_rdc_pos.Y() * 1e3 << std::endl;
    std::cout << "Proton Z pos at RDC: " << proton_rdc_pos.Z() * 1e3 << std::endl;
    std::cout << "Proton aang: " << aang_vec.Theta() / D2R << std::endl;
    std::cout << "E-deg thickenss [mm]: " << edeg_thickness << std::endl;

    csv_data_ofs << tke_initial << "," << proton_theta << ","
                 << proton_phi << "," << "," << proton_rdc_pos.X() * 1e3
                 << "," << proton_rdc_pos.Y() * 1e3
                 << "," << edeg_thickness << std::endl;
  }

  return rdc_xes;
}

std::ostream& s13::sim::operator<<(std::ostream &ostr, const Particle &rhs) {
    ostr << rhs.name() << ": " << std::endl;
    ostr << "\tbeta: " << rhs.beta() << "; gamma: " << rhs.gamma() << std::endl;
    ostr << "\ttke: " << rhs.tke_mev() << " MeV; momentum_mevc: " << rhs.momentum_mevc().Mag() << " MeV" << std::endl;
    //ostr << "\tposition [m]: " << rhs.position().X() << ", " << rhs.position().Y() << ", " << rhs.position().Z();
    ostr << "\tdistance: " << rhs.position().Mag() << " m";
    ostr << "; theta: " << rhs.direction().Theta() * 1/D2R << " deg." << std::endl;

    return ostr;
}


///////////////////////
/// main proc
///////////////////////


int main(int argc, char *argv[]) {
    std::vector<Double_t> eff_angles_range;
    for (Double_t angle = -24; angle <= 24; angle += 0.2) {
        eff_angles_range.push_back(angle);
    }

    // auto measured_angles_range = s13::sim::proton_theta_measured(eff_angles_range, {0,0,-1.40e-3});
    auto measured_xes_range =
      s13::sim::proton_delta_x_vs_phi(eff_angles_range,
                                      55+( 1.47 / 2 ), {0,0,0});
    //std::cout << "effective angle: " << p_theta_eff << std::endl;

    return 0;
}
