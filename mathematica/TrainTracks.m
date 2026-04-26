BeginPackage["TrainTracks`"]

TrainTracks::usage = "Functions for manipulating and drawing train tracks and their automata."

PlotAutomaton::usage = "Plot a train track automaton."

PlotTrainTrack::usage = "Plot a train track."

PlotAutomatonTrainTracks::usage = "Plot train tracks in an automaton as a table."

DropSelfLoops::usage = "Drop edges that connect a vertex to itself in an automaton."

(* Options for PlotAutomaton *)

PlotTrainTracks::usage = "PlotTrainTracks is an option for PlotAutomaton to plot train track at each vertex in automaton [Default True]."

TrainTrackOffset::usage = "TrainTrackOffset is an option for PlotAutomaton to offset each train track from centre of vertex [Default {0,0}]."

TrainTrackSize::usage = "TrainTrackSize is an option for PlotAutomaton to set the size of plotted train tracks [Default 0.25]."


Begin["`Private`"]

PlotTrainTrack[tt_, opts:OptionsPattern[]] :=
  GraphPlot[tt, opts, VertexLabeling -> Automatic, PlotStyle -> Thick,
    VertexRenderingFunction -> ({White, EdgeForm[Black], Disk[#, .1],
                                 Black} &)]

(* PlotTrainTrack inherits the options for GraphPlot *)
Options[PlotTrainTrack] = Options[GraphPlot]


PlotAutomatonTrainTracks[g_, opts:OptionsPattern[]] :=
  PlotTrainTrack[#, opts] & /@ (#[[2]] & /@ g[[1]])

(* PlotAutomatonTrainTrack inherits the options for GraphPlot *)
Options[PlotAutomatonTrainTrack] = Options[GraphPlot]


PlotAutomaton[g_, opts:OptionsPattern[]] :=
  If[OptionValue[PlotTrainTracks] == True,
     (* Draw train track at each vertex with coding as tooltip. *)
     Module[{tt = PlotAutomatonTrainTracks[g]},
       GraphPlotAutomaton[g, FilterRules[{opts},Options[GraphPlotAutomaton]],
         VertexRenderingFunction ->
           (Tooltip[Inset[Framed[tt[[#2]], Background -> LightYellow],
            #1+OptionValue[TrainTrackOffset], {0,0},
            OptionValue[TrainTrackSize]],
            TrainTrackCodingString[g,#2]] &)
       ]
     ],
     (* Print a disk at each vertex with coding as tooltip. *)
     GraphPlotAutomaton[g, FilterRules[{opts},Options[GraphPlotAutomaton]],
       VertexRenderingFunction ->
           (Tooltip[{Brown, EdgeForm[Black], Disk[#, .025],
                    Black}, TrainTrackCodingString[g,#2]] &)
     ]
  ]

(* PlotAutomaton inherits the options for GraphPlot *)
Options[PlotAutomaton] =
  Join[
    {PlotTrainTracks -> False, TrainTrackOffset -> {0,0},
     TrainTrackSize -> .25},
    Options[GraphPlot]
  ]

DropSelfLoops[g_] := {g[[1]], Select[g[[2]], #[[1, 1]] != #[[1, 2]] &]}

(*
   Helper functions (private context)
*)

GraphPlotAutomaton[g_, ex___] :=
  If[Length[g[[2]]] <= 100,
    GraphPlot[#[[1]] & /@ g[[2]], ex, DirectedEdges -> True,
      VertexLabeling -> True, PlotStyle -> Thick, SelfLoopStyle -> .25],
    If[Length[g[[2]]] <= 200,
      GraphPlot[#[[1]] & /@ g[[2]], ex, DirectedEdges -> False,
        VertexLabeling -> False, PlotStyle -> Thick, SelfLoopStyle -> .1],
      GraphPlot[#[[1]] & /@ g[[2]], ex, DirectedEdges -> False,
        VertexLabeling -> False, SelfLoopStyle -> .1]
    ]
 ]

(* GraphPlotAutomaton inherits the options for GraphPlot *)
Options[GraphPlotAutomaton] = Options[GraphPlot]

(* A string with the train track # and its coding (for tooltips) *)
TrainTrackCodingString[g_, idx_] :=
  StringJoin[ToString[idx],": ", g[[1, idx, 1]]]

End[(* "`Private`" *)]

EndPackage[(* "TrainTracks`" *)]
