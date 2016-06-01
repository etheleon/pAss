#!/usr/bin/env perl

use strict;
use warnings;

ok `which xvfb-run`, "xvfb is installed";
ok `which muscle`, "Muscle is installed";
ok `which MEGAN`, "MEGAN is installed";
