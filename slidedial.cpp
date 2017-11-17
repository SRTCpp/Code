#include "slidedial.h"

Jcrap::slidedial::slidedial(string labelstr, QWidget *parent, const char *name)
	: QHBox(parent, name), 
	labeller (labelstr.c_str(), this),
	spinner(this), 
	slider(Qt::Horizontal, this)
{
	setMargin(6);
	setSpacing(6);
	
	rerange(0, 0);
//	setMinimumSize(200, 30);

	QObject::connect(&spinner, SIGNAL(valueChanged(int)),
		&slider, SLOT(setValue(int)));
	QObject::connect(&slider, SIGNAL(valueChanged(int)),
		&spinner, SLOT(setValue(int)));
	QObject::connect(&slider, SIGNAL(valueChanged(int)),
		this, SIGNAL(valueChanged(int)));
	spinner.setValue(0);	
}

void Jcrap::slidedial::value(int i)
{
	spinner.setValue(i);
	
}

void Jcrap::slidedial::rerange(int i, int j)
{
	spinner.setRange(i,j);
	slider.setRange(i,j);
}

int Jcrap::slidedial::value()
{
	return slider.value();
}
