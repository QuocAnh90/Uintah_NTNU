catch {rename DaveW_Writers_SegFldWriter ""}

itcl_class DaveW_Writers_SegFldWriter {
    inherit Module
    method modname {} {
	set n $this
	if {[string first "::" "$n"] == 0} {
	    set n "[string range $n 2 end]"
	}
	return $n
    }
    constructor {config} {
	set name SegFldWriter
	set_defaults
    }
    method set_defaults {} {
	global $this-filetype
	set $this-filetype Binary
    }
    method ui {} {
	set w .ui[modname]
	if {[winfo exists $w]} {
	    raise $w
	    return;
	}
	toplevel $w

	make_labeled_radio $w.filetype "Format:" "" left $this-filetype \
		{Binary ASCII}
	pack $w.filetype
	entry $w.f -textvariable $this-filename -width 40 \
		-borderwidth 2 -relief sunken
	pack $w.f -side bottom
	bind $w.f <Return> "$this-c needexecute "
    }
}
