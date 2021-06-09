library(DiagrammeR)

grViz("digraph flowchart {
      # node definitions with substituted label text
      node [fontname = Helvetica, shape = rectangle]        
      tab1 [label = '@@1']
      tab2 [label = '@@2']
      tab3 [label = '@@3']
      tab4 [label = '@@4']
      tab5 [label = '@@5']

      # edge definitions with the node IDs
      tab1 -> tab2 -> tab3 -> tab4 -> tab5;
      }

      [1]: 'Questionnaire sent to n=1000 participants'
      [2]: 'Participants responded to questionnaire n=850'
      [3]: 'Participants came to clinic for evaluation n=700'
      [4]: 'Participants eligible for the study n=600'
      [5]: 'Study sample n=600'
")


grViz("
digraph {

  # graph attributes
  graph [overlap = true]

  # node attributes
  node [shape = box,
        fontname = Helvetica,
        color = blue]

  # edge attributes
  edge [color = gray]

  # node statements
  A; B; C; D; E
  F [color = black]

  # node attributes
  node [shape = circle,
        fixedsize = true,
        width = 0.9]

  # node statements
  1; 2; 3; 4; 5; 6; 7; 8

  # edge statements
  A->1; B->2                   // gray
  B->3 [color = red]           // red
  B->4                         // gray
  C->A [color = green]         // green
  1->D; E->A; 2->4; 1->5; 1->F // gray
  E->6; 4->6; 5->7; 6->7       // gray
  3->8 [color = blue]          // blue
}
")


grViz("
digraph {

  # graph attributes
  graph [overlap = true]


  
  # node attributes
  node [shape = oval,
        width = 2.5]
  # node statements
  A [label = 'Demographic history']
  Recent [label = 'Recent history']
  Deep [label = 'Deeper in past history']
  MS [label = 'Demographic simulation']
  
  # node methods completed
  node [shape = box,
        fontname = Helvetica,
        style = filled,
        fillcolor = LavenderBlush]
  Lai [label = 'Local ancestry\nXGMix']
  
  # node methods working on
  node [shape = box,
        fontname = Helvetica,
        color = DimGray,
        style = filled,
        fillcolor = Lavender]
  Lai [label = 'Local ancestry\nXGMix']
  Tracts [label = 'Admixture dynamics\n(Tracts)']
  SFS [label = 'Compute joint SFS']
  
  # node methods future
  node [shape = box,
        fontname = Helvetica,
        color = DimGray,
        style = filled,
        fillcolor = Azure]
  joinSFS [label = 'Joint demographic history\n(moments / dadi)']
  integrateDemo [label = 'Integrate inferences']
  simMeztizo [label = 'Simulate mexican admixed\npopulation (msprime)']
  MLS [label = 'Simuation of mutation load\n(SLiM, add selection)']
  RA [label = 'Daleterious variants\nrecesive or additive?']
  # edge attributes
  edge [color = gray]

  # edge statements
  A->Recent                   // gray
  A->Deep                    // gray
  Recent->Lai->Tracts
  Deep->SFS->joinSFS
  Tracts->integrateDemo
  joinSFS->integrateDemo
  integrateDemo->MS->simMeztizo
  MS->MLS
}

")

