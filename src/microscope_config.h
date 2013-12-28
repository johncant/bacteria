#ifndef __BACTERIA__MICROSCOPE_CONFIG_H__
#define __BACTERIA__MICROSCOPE_CONFIG_H__

namespace Bacteria {

  class MicroscopeConfig {
    public:
      double scale;
      double fg_intensity, bg_intensity;
      int view_width, view_height;
  };

}

#endif
