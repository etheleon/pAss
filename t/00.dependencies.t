#!/usr/bin/env perl

use strict;
use warnings;
use Test::More;
use IPC::Cmd qw/can_run run/;


ok (can_run("xvfb-run"), "xvfb is installed");
ok (can_run("muscle"), "Muscle is installed");
ok (can_run("MEGAN"), "MEGAN is installed");

done_testing;
