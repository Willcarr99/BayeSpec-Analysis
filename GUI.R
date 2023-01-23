## The GUI (tcltk)
##---------------------------------------------------------------

library("tcltk")

                                        # Extra functions
GUIopenfile <- function(){
  fileName <- tclvalue(tkgetOpenFile(filetypes=
                                     "{{ascii spectrum} {.dat .asc}}"))
  load.dat(fileName)
}

GUIopenfilecomp <- function(){
  fileName <- tclvalue(tkgetOpenFile(filetypes=
                                     "{{ascii spectrum} {.dat .asc}}"))
  overlay(fileName)
}


##############################################
## The main window
GUIwindow <- function(){
  tt <<- tktoplevel()

  tktitle(tt) <- "BaySpec Commands"

                                        # Buttons
  b.open <- tkbutton(tt, text="Open",command=GUIopenfile)

  b.print <- tkbutton(tt, text="Print",command=print.spec)

  b.loadcal <- tkbutton(tt, text="Load Calibration",command=load.cal)
  b.savecal <- tkbutton(tt, text="Save Calibration",command=savecal)

  b.lines <- tkbutton(tt, text="Load Lines",command=load.lines)
  
  b.log <- tkbutton(tt, text="Log/Lin",command=logy)

  b.goto <- tkbutton(tt, text="Goto",command=goto)

  b.get <- tkbutton(tt, text="Get Point",command=getpoint)

  b.zoomx <- tkbutton(tt, text="Zoom X",command=zx)

  b.zoomox <- tkbutton(tt, text="Zoom out",command=zo)

  b.zoomy <- tkbutton(tt, text="Zoom Y",command=zy)

  b.zoomay <- tkbutton(tt, text="Auto Zoom Y",command=ay)

  b.zoomall <- tkbutton(tt, text="Auto Zoom All",command=za)

  b.gross <- tkbutton(tt, text="Gross Area",command=gross)

  b.net <- tkbutton(tt, text="Net Area",command=net)

  b.ul <- tkbutton(tt, text="Upper Limit",command=ul)

  b.cal <- tkbutton(tt, text="Calibrate",command=cal)

  b.rebin <- tkbutton(tt, text="Rebin",command=rebin)

  b.refresh <- tkbutton(tt, text="Refresh",command=replot)

  b.invert <- tkbutton(tt, text="Invert",command=invert)

  b.compare <- tkbutton(tt, text="Overlay",command=GUIopenfilecomp)
  
  b.quit <- tkbutton(tt, text="Quit",command=quit)

                                        # Pack everything together and display
                                        #tkpack(b.open, b.print, b.loadcal, b.savecal, b.log, b.goto, b.zoomx,
                                        #       b.zoomox,b.zoomy, b.zoomay, b.zoomall,
                                        #       b.gross, b.net, b.ul, b.cal, b.rebin, b.refresh, b.invert, b.quit)
  i <- 1
  tkgrid(b.open,row=i,column=1)
  i <- i+1
  tkgrid(b.print,row=i,column=1)
  i <- i+1
  tkgrid(b.loadcal,row=i,column=1)
  i <- i+1
  tkgrid(b.savecal,row=i,column=1)
  i <- i+1
  tkgrid(b.cal,row=i,column=1)
  i <- i+1
  tkgrid(b.lines,row=i,column=1)
  i <- i+1
  tkgrid(b.gross,row=i,column=1)
  i <- i+1
  tkgrid(b.net,row=i,column=1)
  i <- i+1
  tkgrid(b.ul,row=i,column=1)
  i <- i+1
  tkgrid(b.compare,row=i,column=1)
  i <- i+1
  tkgrid(b.quit,row=i,column=1)

  i <- 1
  tkgrid(b.log,row=i,column=2)
  i <- i+1
  tkgrid(b.goto,row=i,column=2)
  i <- i+1
  tkgrid(b.get,row=i,column=2)
  i <- i+1
  tkgrid(b.zoomx,row=i,column=2)
  i <- i+1
  tkgrid(b.zoomy,row=i,column=2)
  i <- i+1
  tkgrid(b.zoomay,row=i,column=2)
  i <- i+1
  tkgrid(b.zoomox,row=i,column=2)
  i <- i+1
  tkgrid(b.zoomall,row=i,column=2)
  i <- i+1
  tkgrid(b.rebin,row=i,column=2)
  i <- i+1
  tkgrid(b.invert,row=i,column=2)
  i <- i+1
  tkgrid(b.refresh,row=i,column=2)

  tkwm.geometry(tt,"-5+5")

}


GUIwindow()
