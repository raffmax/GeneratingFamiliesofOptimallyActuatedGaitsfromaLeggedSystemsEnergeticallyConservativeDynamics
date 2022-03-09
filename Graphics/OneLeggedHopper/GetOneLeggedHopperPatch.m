% *************************************************************************
%
% function [f, v, c] = GetOneLeggedHopperPatch(x, y, phi, alpha, l, ualpha, ul)
% 
% Returns faces, vertices, and color-values representing a patch that shows
% the 3D view of a monopod hopper with two degrees of freedom.
%
% GetMonopod2DOFPatch(x, y, phi, alpha, ualpha, l, ul), without return
% arguments shows the patch in a new figure
%
% Input:  - x, y, phi (position and orientation of the main body)
%         - alpha, l (length and orientation of the leg)
%         - ualpha, ul (deflection of the actuators from their resting
%           position)
%
% Output: - faces 'f', vertices 'v', and color values 'c' as used in the
%           'patch' function
%
% Created by C. David Remy on 03/14/2011
% MATLAB 2010a
%
% Documentation:
%  'A MATLAB Framework For Gait Creation', 2011, C. David Remy (1), Keith
%  Buffinton (2), and Roland Siegwart (1),  International Conference on
%  Intelligent Robots and Systems, September 25-30, San Francisco, USA 
%
% (1) Autonomous Systems Lab, Institute of Robotics and Intelligent Systems, 
%     Swiss Federal Institute of Technology (ETHZ) 
%     Tannenstr. 3 / CLA-E-32.1
%     8092 Zurich, Switzerland  
%     cremy@ethz.ch; rsiegwart@ethz.ch
%
% (2) Department of Mechanical Engineering, 
%     Bucknell University
%     701 Moore Avenue
%     Lewisburg, PA-17837, USA
%     buffintk@bucknell.edu
%
%   See also OUTPUTCLASS, PATCH.
%
function [f, v, c] = GetOneLeggedHopperPatch(x,y,phi,alpha,l,epsilon,edges,colors)

    % The faces, vertices, and color values are only created once and then
    % stored in memory:
    persistent mainBodyVertices upperLegVertices lowerLegVertices upperSpringVertices lowerSpringVertices upperActuatorVertices lowerActuatorVertices
    persistent f1 f2 f3 f4 f5 f6 f7
    persistent c1 c2 c3 c4 c5 c6 c7
    %if isempty(mainBodyVertices) || isempty(upperLegVertices) || isempty(lowerLegVertices) || isempty(upperSpringVertices) || isempty(lowerSpringVertices) || isempty(upperActuatorVertices) || isempty(lowerActuatorVertices)
        [f1, mainBodyVertices, c1]      = GetMainBodyMonopodPatch(colors,edges,epsilon);
        [f2, upperLegVertices, c2]      = GetUpperLegPatch(colors,epsilon);
        [f3, lowerLegVertices, c3]      = GetLowerLegPatch(colors,epsilon);
        [f4, upperSpringVertices, c4]   = GetUpperSpringPatch(colors.spring);
        [f5, lowerSpringVertices, c5]   = GetLowerSpringPatch(colors.spring);
        [f7, upperActuatorVertices, c7] = GetUpperActuatorPatch(colors);
        [f6, lowerActuatorVertices, c6] = GetLowerActuatorPatch(colors);
    %end
    
    % NOTE: Everything is implicitly transformed from the [X Y Z] in which
    % the dynamics are stated into the [X -Z Y] system in which the
    % graphics are displayed, i.e. the translation of the main body is
    % [x,0,y] instead of [x,y,0]. The main body and leg angles are thus
    % inverted: 
    phi    = -phi;
    alpha  = -alpha;
           
    %% Transform the vertices, and compose the complete patch object
    %% Main body:
    v = TransformVertices(mainBodyVertices, Body312dc([0,0,phi]), [x,0,y]);f = f1;c = c1;
    %% Upper leg segments (rotate and shift side wards)
    v_ = TransformVertices(upperLegVertices,Body312dc([0,0,alpha]),[0,0,0]); % rotate
    v2 = TransformVertices(v_,diag([1,+1,1])*Body312dc([0,0,phi]), [x,+0.6,y]); % shift
    [v,f,c]=AddPatchesWithColor(v,f,c,v2,f2,c2);
    v2 = TransformVertices(v_,diag([1,-1,1])*Body312dc([0,0,phi]), [x,-0.6,y]); % shift
    [v,f,c]=AddPatchesWithColor(v,f,c,v2,f2,c2);
    %% Lower leg segments (move downwards, rotate and shift side wards)
    v_ = TransformVertices(lowerLegVertices,eye(3),[0,0,-l-0.05*(1-1)]);
    v_ = TransformVertices(v_,Body312dc([0,0,alpha]),[0,0,0]); % rotate
    v3 = TransformVertices(v_,Body312dc([0,0,phi]), [x,+0.6,y]); % shift 
    [v,f,c]=AddPatchesWithColor(v,f,c,v3,f3,c3);
    v3 = TransformVertices(v_,Body312dc([0,0,phi]), [x,-0.6,y]); % shift
    [v,f,c]=AddPatchesWithColor(v,f,c,v3,f3,c3);
    %% Lower actuators (move downwards, rotate and shift side wards)
    v_ = TransformVertices(lowerActuatorVertices,eye(3),[0,0,(-0.5+0.05)]);
    v_ = TransformVertices(v_,Body312dc([0,0,alpha]),[0,0,0]); % rotate
    v6 = TransformVertices(v_,Body312dc([0,0,phi]), [x,+0.6,y]); % shift 
    [v,f,c]=AddPatchesWithColor(v,f,c,v6,f6,c6);
    v6 = TransformVertices(v_,Body312dc([0,0,phi]), [x,-0.6,y]); % shift  
    [v,f,c]=AddPatchesWithColor(v,f,c,v6,f6,c6);              
    %% Lower spring (scale, move downwards, rotate and shift side wards)
    delta_l = (1-(0.95-l-0.05*(1-1))/0.34);
    v_ = TransformVertices(lowerSpringVertices,diag([1,1,delta_l]),[0,0,(-l + 0.16-0.05*(1-1))]);
    v_ = TransformVertices(v_,Body312dc([0,0,alpha]),[0,0,0]); % rotate
    v5 = TransformVertices(v_,Body312dc([0,0,phi]), [x,+0.6,y]); % shift 
    [v,f,c]=AddPatchesWithColor(v,f,c,v5,f5,c5);
    v5 = TransformVertices(v_,Body312dc([0,0,phi]), [x,-0.6,y]); % shift 
    [v,f,c]=AddPatchesWithColor(v,f,c,v5,f5,c5);
    %% Upper actuators (rotate and shift side wards)
    v_ = TransformVertices(upperActuatorVertices,Body312dc([0,0,0]),[0,0,0]); % rotate
    v7 = TransformVertices(v_,Body312dc([0,0,phi]), [x,+0.6,y]); % shift 
    [v,f,c]=AddPatchesWithColor(v,f,c,v7,f7,c7);
    v7 = TransformVertices(v_,Body312dc([0,0,phi]), [x,-0.6,y]); % shift  
    [v,f,c]=AddPatchesWithColor(v,f,c,v7,f7,c7);              
    %% Upper springs (bent& compress, rotate and shift side wards)
    delta_alpha = 1.57 -alpha;
    % Bend including compression:
    v_ = BendVertices(upperSpringVertices,0.25/(delta_alpha)*2*pi);
    v4 = TransformVertices(v_,Body312dc([0,0,phi]), [x,+0.6,y]); % shift
    [v,f,c]=AddPatchesWithColor(v,f,c,v4,f4,c4);
    v4 = TransformVertices(v_,Body312dc([0,0,phi]), [x,-0.6,y]); % shift 
    [v,f,c]=AddPatchesWithColor(v,f,c,v4,f4,c4);
    
    % If no output argument is required, this function simply puts the
    % given configuration in a figure:
    if nargout == 0
        figure;
        p = patch('faces', f, 'vertices' ,v,'FaceVertexCData',c,'FaceColor','flat');
        view(3);
        axis equal;
        set(p, 'FaceLighting','phong');
        set(p, 'FaceColor','flat');
        set(p, 'EdgeColor','none'); 
        camlight right
    end
end
% *************************************************************************
% *************************************************************************