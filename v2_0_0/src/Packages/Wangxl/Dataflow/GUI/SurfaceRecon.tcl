itcl_class Wangxl_Surface_SurfaceRecon {
    inherit Module
    constructor {config} {
        set name SurfaceRecon
        set_defaults
    }

    method set_defaults {} {
    }

    method ui {} {
        set w .ui[modname]
        if {[winfo exists $w]} {
            raise $w
            return
        }
        toplevel $w
        label $w.row1 -text "This GUI was auto-generated by the Component Wizard."
        label $w.row2 -text {edit the file "/home/sci/wangxl/SCIRun/src/Packages/Wangxl/Dataflow/GUI/SurfaceRecon.tcl" to modify it.}
        pack $w.row1 $w.row2 -side top -padx 10 -pady 10
    }
}


