{

    /*
      auto ds = s13::ana::MakeDataset({273});
      TTree *tt = ds.GetTree();
    */

    branch_list = tt->GetListOfBranches();
    for (auto *entry : *branch_list) {
        //cout << "entry no.: " << i << " pointer: " << entry << endl;
        cout << entry->GetName() << endl;
    }
}
