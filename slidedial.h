#include <qslider.h>
#include <qlabel.h>
#include <qlcdnumber.h>
#include <qspinbox.h>
#include <q3hbox.h>
#include <q3vbox.h>
#include <q3grid.h>
#include <qlineedit.h>
#include <q3buttongroup.h>
#include <qradiobutton.h>

#include "Jcrap.h"

#ifndef SLIDEDIAL_H
#define SLIDEDIAL_H 1

class Jcrap::slidedial : public Q3HBox
{
	Q_OBJECT

public:
	slidedial(string labelstr="A Slidedial", QWidget *parent=0, const char *name=0);
	void rerange(int, int);
	int value();
	void value(int);
	int upperlimit();

signals:
	void valueChanged( int );

private:
	QLabel   labeller;
	QSpinBox spinner;
	QSlider  slider;
};

#endif
