/**
 *   \file rootscript.hpp
 *   \brief Declaration of RootScript - convenience object for drawing outputs...
 *
 */

#pragma once

#include <iostream>

#include <TStyle.h>
#include <TString.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TFile.h>


namespace s13 {
  namespace ana {

    enum class PadSequence { row, column };
    using PadSeq = enum PadSequence;


    class RootScript {

    public:

      explicit RootScript(TString name, TString out_dir_name = "out") :
        m_name{name}, m_is_first_page{true},
        m_is_page_printed{true}, m_canvas_idx{1},
        m_pad_seq{PadSeq::column} {
          Prepare(out_dir_name);
        }

      ~RootScript() {
        Close();
      }

      TString GetName() const {
        return m_name;
      }

      void SetName(TString name) {
        m_name = name;
      }

      TCanvas& GetCanvas() {
        return *m_root_canvas;
      }

      void ActivateCanvas() {
        m_root_canvas->cd();
      }

      TChain* GetChain() const {
        return m_tchain;
      }

      std::vector<TChain*> GetChainSegments() {
        if (m_tchains_vec.empty()) {
          return std::vector<TChain*>{m_tchain};
        } else {
          return m_tchains_vec;
        }
      }

      void cd(PadSeq opt) {
        m_pad_seq = opt;
      }

      void cd(int pad_id = -1) {
        if (pad_id == -1) {
          m_canvas_idx += 1;
          //std::cout << "m_canvas_idx = " << m_canvas_idx << std::endl;
          if (m_pad_seq == PadSeq::row) {
            m_root_canvas->cd(m_canvas_idx);
          } else {
            int col = m_canvas_idx / m_page_rows;
            int row = m_canvas_idx % m_page_rows;
            if (row == 0) {
              row = m_page_rows;
            } else {
              col += 1;
            }
            int idx = (row - 1) * m_page_columns + col;
            // std::cout << "setting next pad to col. " << col << " row "
            //           << row << " with idx: " << idx << std::endl;
            m_root_canvas->cd(idx);
          }
        } else {
          m_root_canvas->cd(pad_id);
          m_canvas_idx = pad_id;
        }
      }

      RootScript& NewPage(int col_num = 1, int row_num = 1) {
        m_page_rows = row_num;
        m_page_columns = col_num;
        PrintPage();
        m_canvas_idx = 0;
        m_root_canvas->Clear();
        m_root_canvas->Divide(col_num, row_num);
        m_root_canvas->cd(0);
        m_is_page_printed = false;
        return *this;
      }

      void Prepare(TString out_dir_name = "out") {
        m_out_dir_name = out_dir_name;
        m_root_filename = TString::Format("%s/%s.root",
                                          m_out_dir_name.Data(),
                                          m_name.Data());
        m_pdf_filename = TString::Format("%s/%s.pdf",
                                         m_out_dir_name.Data(),
                                         m_name.Data());
        // m_root_file = std::unique_ptr<TFile>(new TFile(m_root_filename.Data(), "RECREATE"));
        m_root_canvas = std::unique_ptr<TCanvas>(new TCanvas(m_name.Data(), m_name.Data()));
        gStyle->SetOptStat(1111111);
      }

      void Close() {
        // close pdf output file
        if (m_is_first_page) {
          m_root_canvas->Print(m_pdf_filename, "pdf");
        } else {
          m_root_canvas->Print(m_pdf_filename + ")", "pdf");
        }

        // KLUDGE: .root files were rarely used, stop support its writing
        // write & close root file
        // m_root_file->Write();
        // m_root_file->Close();
      }

      TString GetFilename(TString suffix = "") {
        return s13::ana::tstrfmt("%s/%s%s", m_out_dir_name.Data(),
                                 m_name.Data(), suffix.Data());
      }

    private:

      void PrintPage() {
        if (m_is_page_printed) {
          return;
        }

        if (m_is_first_page) {
          m_root_canvas->Print(m_pdf_filename + "(", "pdf");
        } else {
          m_root_canvas->Print(m_pdf_filename, "pdf");
        }

        m_is_page_printed = true;
        m_is_first_page = false;
      }

      TString m_name;
      TString m_out_dir_name;
      TString m_root_filename;
      TString m_pdf_filename;
      bool m_is_first_page;
      bool m_is_page_printed;
      int m_canvas_idx;
      int m_page_rows;
      int m_page_columns;
      PadSeq m_pad_seq;
      std::unique_ptr<TFile> m_root_file;
      std::unique_ptr<TCanvas> m_root_canvas;
      TChain* m_tchain;
      std::vector<TChain*> m_tchains_vec;
    };

  }
}
