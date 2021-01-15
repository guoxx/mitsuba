#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/hw/basicshader.h>
#include <mitsuba/core/warp.h>

//#pragma optimize("", off)

MTS_NAMESPACE_BEGIN

class Disney : public BSDF {

    struct SurfaceData
    {
        Spectrum baseColor;
        Spectrum specularColor;
        float subsurface;
        float metallic;
        float specular;
        float specularTint;
        float roughness;
        float anisotropic;
        float sheen;
        float sheenTint;
        float clearcoat;
        float clearcoatGloss;
    };

public:
    Disney(const Properties &props)
        : BSDF(props) {

        m_baseColor = new ConstantSpectrumTexture(props.getSpectrum("baseColor", Spectrum(.5f)));
        m_subsurface = new ConstantFloatTexture(props.getFloat("subsurface", .0f));
        m_metallic = new ConstantFloatTexture(props.getFloat("metallic", .0f));
        m_specular = new ConstantFloatTexture(props.getFloat("specular", .5f));
        m_specularTint = new ConstantFloatTexture(props.getFloat("specularTint", .0f));
        m_roughness = new ConstantFloatTexture(props.getFloat("roughness", .5f));
        m_anisotropic = new ConstantFloatTexture(props.getFloat("anisotropic", .0f));
        m_sheen = new ConstantFloatTexture(props.getFloat("sheen", .0f));
        m_sheenTint = new ConstantFloatTexture(props.getFloat("sheenTint", .5f));
        m_clearcoat = new ConstantFloatTexture(props.getFloat("clearcoat", .0f));
        m_clearcoatGloss = new ConstantFloatTexture(props.getFloat("clearcoatGloss", 1.f));
    }

    Disney(Stream *stream, InstanceManager *manager)
        : BSDF(stream, manager) {
        m_baseColor = static_cast<Texture *>(manager->getInstance(stream));
        m_subsurface = static_cast<Texture *>(manager->getInstance(stream));
        m_metallic = static_cast<Texture *>(manager->getInstance(stream));
        m_specular = static_cast<Texture *>(manager->getInstance(stream));
        m_specularTint = static_cast<Texture *>(manager->getInstance(stream));
        m_roughness = static_cast<Texture *>(manager->getInstance(stream));
        m_anisotropic = static_cast<Texture *>(manager->getInstance(stream));
        m_sheen = static_cast<Texture *>(manager->getInstance(stream));
        m_sheenTint = static_cast<Texture *>(manager->getInstance(stream));
        m_clearcoat = static_cast<Texture *>(manager->getInstance(stream));
        m_clearcoatGloss = static_cast<Texture *>(manager->getInstance(stream));

        configure();
    }

    void serialize(Stream *stream, InstanceManager *manager) const override {
        BSDF::serialize(stream, manager);

        manager->serialize(stream, m_baseColor.get());
        manager->serialize(stream, m_subsurface.get());
        manager->serialize(stream, m_metallic.get());
        manager->serialize(stream, m_specular.get());
        manager->serialize(stream, m_specularTint.get());
        manager->serialize(stream, m_roughness.get());
        manager->serialize(stream, m_anisotropic.get());
        manager->serialize(stream, m_sheen.get());
        manager->serialize(stream, m_sheenTint.get());
        manager->serialize(stream, m_clearcoat.get());
        manager->serialize(stream, m_clearcoatGloss.get());
    }

    void configure() override {
        m_components.clear();
        m_components.push_back(EGlossyReflection | EAnisotropic | EFrontSide | ESpatiallyVarying);
        m_components.push_back(EDiffuseReflection | EAnisotropic | EFrontSide | ESpatiallyVarying);

        m_usesRayDifferentials = m_baseColor->usesRayDifferentials() ||
                                 m_subsurface->usesRayDifferentials() ||
                                 m_metallic->usesRayDifferentials() ||
                                 m_specular->usesRayDifferentials() ||
                                 m_specularTint->usesRayDifferentials() ||
                                 m_roughness->usesRayDifferentials() ||
                                 m_anisotropic->usesRayDifferentials() ||
                                 m_sheen->usesRayDifferentials() ||
                                 m_sheenTint->usesRayDifferentials() ||
                                 m_clearcoat->usesRayDifferentials() ||
                                 m_clearcoatGloss->usesRayDifferentials();

        BSDF::configure();
    }

    Spectrum getDiffuseReflectance(const Intersection &its) const override {
        NotImplementedError("getDiffuseReflectance");
        return Spectrum(0.0f);
    }

    Float sqr(Float x) const { return x * x; }

    Float roughnessToAlpha(Float roughness) const { return roughness * roughness; }

    // https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
    //
    // The Schlick Fresnel approximation is:
    //
    // R = R(0) + (1 - R(0)) (1 - cos theta)^5,
    //
    // where R(0) is the reflectance at normal indicence.
    Float schlickWeight(Float cosTheta) const {
        auto m = math::clamp(1.0f - cosTheta, 0.0f, 1.0f);
        return (m * m) * (m * m) * m;
    }

    Float FrSchlick(Float R0, Float cosTheta) const {
        return math::lerp(schlickWeight(cosTheta), R0, 1.0f);
    }

    Spectrum FrSchlick(const Spectrum &R0, Float cosTheta) const {
        return math::lerp(Spectrum(schlickWeight(cosTheta)), R0, Spectrum(1.));
    }

    Float GTR1(Float NdotH, Float a) const
    {
        if (a >= 1)
            return INV_PI;

        auto a2 = a * a;
        auto t = 1 + (a2 - 1) * NdotH * NdotH;
        return (a2 - 1) / (M_PI * std::log(a2) * t);
    }

    Vector3f sampleGTR1(const Point2f& sample_, Float alpha) const
    {
        auto a2 = alpha * alpha;
        auto cosTheta2 = (1 - std::pow(a2, 1 - sample_.x)) / (1 - a2);
        auto cosTheta = math::safe_sqrt(cosTheta2);
        auto sinTheta = math::safe_sqrt(1 - cosTheta2);
        auto phi = sample_.y * 2 * M_PI;
        return Vector3f(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta);
    }

    Float pdfGTR1(const Vector3f& v, Float alpha) const
    {
        if (Frame::cosTheta(v) > 0)
        {
            const auto cosTheta = Frame::cosTheta(v);
            return GTR1(cosTheta, alpha) * cosTheta;
        }
        return 0.0f;
    }

    void stretchAniso(Float anisotropic, Float roughness, Float* alphaU, Float* alphaV) const
    {
        auto aspect = math::safe_sqrt(1 - anisotropic * .9f);
        auto ax = std::max(.001f, roughnessToAlpha(roughness) / aspect);
        auto ay = std::max(.001f, roughnessToAlpha(roughness) * aspect);        
        *alphaU = ax;
        *alphaV = ay;
    }

    Float GTR2_aniso(const Vector3f& v, Float anisotropic, Float roughness) const
    {
        auto cosTheta = v.z;
        auto sinThetaXcosPhi = v.x;
        auto sinThetaXsinPhi = v.y;

        Float ax, ay;
        stretchAniso(anisotropic, roughness, &ax, &ay);

        return 1 / (M_PI * ax * ay * sqr(sqr(sinThetaXcosPhi / ax) + sqr(sinThetaXsinPhi / ay) + sqr(cosTheta)));
    }

    Vector3f sampleGTR2Aniso(const Point2f& sample_, Float anisotropic, Float roughness) const
    {
        Float ax, ay;
        stretchAniso(anisotropic, roughness, &ax, &ay);

        Vector3f vx = Vector3f(1, 0, 0);
        Vector3f vy = Vector3f(0, 1, 0);
        Vector3f n = Vector3f(0, 0, 1);
        auto f = math::safe_sqrt(sample_.y / (1 - sample_.y));
        Vector3f x = f * ax * std::cos(2 * M_PI * sample_.x) * vx;
        Vector3f y = f * ay * std::sin(2 * M_PI * sample_.x) * vy;
        Vector3f hp = x + y + n;
        return normalize(hp);
    }

    Float pdfGTR2Aniso(const Vector3f& v, Float anisotropic, Float roughness) const
    {
        if (Frame::cosTheta(v) > 0.0f)
        {
            const auto cosTheta = Frame::cosTheta(v);
            return GTR2_aniso(v, anisotropic, roughness) * cosTheta;
        }
        return 0.0f;
    }

    Float smithG_GGX(Float NdotV, Float alphaG) const
    {
        auto a = alphaG * alphaG;
        auto b = NdotV * NdotV;
        return 1 / (NdotV + math::safe_sqrt(a + b - a * b));
    }

    Float smithG_GGX_aniso(Float NdotV, Float VdotX, Float VdotY, Float ax, Float ay) const
    {
        return 1 / (NdotV + math::safe_sqrt(sqr(VdotX*ax) + sqr(VdotY*ay) + sqr(NdotV)));
    }

    Spectrum diffuseLobe(const BSDFSamplingRecord &bRec, Float roughness) const
    {
        const Vector wh = normalize((bRec.wi + bRec.wo));
        const auto NdotL = Frame::cosTheta(bRec.wi);
        const auto NdotV = Frame::cosTheta(bRec.wo);
        const auto LdotH = dot(bRec.wi, wh);

        const auto FL = schlickWeight(NdotL);
        const auto FV = schlickWeight(NdotV);

        // Diffuse fresnel - go from 1 at normal incidence to .5 at grazing
        // and mix in diffuse retro-reflection based on roughness
        auto Fd90 = 0.5f + 2.0f * LdotH * LdotH * roughness;
        auto Fd = INV_PI * math::lerp(FL, 1.0f, Fd90) * math::lerp(FV, 1.0f, Fd90);
        return Spectrum(Fd);
    }

    Spectrum sssLobe(const BSDFSamplingRecord &bRec, Float roughness) const
    {
        const Vector wh = normalize((bRec.wi + bRec.wo));
        const auto NdotL = Frame::cosTheta(bRec.wi);
        const auto NdotV = Frame::cosTheta(bRec.wo);
        const auto LdotH = dot(bRec.wi, wh);

        const auto FL = schlickWeight(NdotL);
        const auto FV = schlickWeight(NdotV);

        // Based on Hanrahan-Krueger brdf approximation of isotropic bssrdf
        // 1.25 scale is used to (roughly) preserve albedo
        // Fss90 used to "flatten" retroreflection based on roughness
        auto Fss90 = LdotH * LdotH * roughness;
        auto Fss = math::lerp(FL, 1.0f, Fss90) * math::lerp(FV, 1.0f, Fss90);
        auto ss = INV_PI * 1.25f * (Fss * (1.0f / (NdotL + NdotV) - .5f) + .5f);
        return Spectrum(ss);
    }

    Spectrum specularLobe(const BSDFSamplingRecord &bRec, Float roughness, Float anisotropic, Spectrum specularColor) const
    {
        const Vector wh = normalize((bRec.wi + bRec.wo));
        const auto NdotL = Frame::cosTheta(bRec.wi);
        const auto NdotV = Frame::cosTheta(bRec.wo);
        const auto LdotH = dot(bRec.wi, wh);

        const auto LdotX = bRec.wi.x;
        const auto LdotY = bRec.wi.y;
        const auto VdotX = bRec.wo.x;
        const auto VdotY = bRec.wo.y;

        Float ax, ay;
        stretchAniso(anisotropic, roughness, &ax, &ay);
        auto Ds = GTR2_aniso(wh, anisotropic, roughness);
        Spectrum Fs = FrSchlick(specularColor, LdotH);
        auto Gs = smithG_GGX_aniso(NdotL, LdotX, LdotY, ax, ay);
        Gs *= smithG_GGX_aniso(NdotV, VdotX, VdotY, ax, ay);

        return Ds * Fs * Gs;
    }

    Spectrum sheenLobe(const BSDFSamplingRecord &bRec, Float sheen, Spectrum sheenColor) const
    {
        const Vector wh = normalize((bRec.wi + bRec.wo));
        const auto LdotH = dot(bRec.wi, wh);

        const auto FH = schlickWeight(LdotH);
        return FH * sheen * sheenColor;
    }

    Spectrum clearcoatLobe(const BSDFSamplingRecord &bRec, Float clearcoat, Float clearcoatGloss) const
    {
        const Vector wh = normalize((bRec.wi + bRec.wo));
        const auto NdotL = Frame::cosTheta(bRec.wi);
        const auto NdotV = Frame::cosTheta(bRec.wo);
        const auto LdotH = dot(bRec.wi, wh);
        const auto NdotH = Frame::cosTheta(wh);

        // clearcoat (ior = 1.5 -> F0 = 0.04)
        auto Dr = GTR1(NdotH, math::lerp<Float>(clearcoatGloss, .1f, .001f));
        auto Fr = FrSchlick(.04f, LdotH);
        auto Gr = smithG_GGX(NdotL, .25) * smithG_GGX(NdotV, .25);
        return Spectrum(0.25f * clearcoat * Gr * Fr * Dr);
    }

    SurfaceData prepareSurfaceData(const Intersection& its) const {
        SurfaceData sd;
        sd.baseColor = m_baseColor->eval(its);
        sd.subsurface = m_subsurface->eval(its).average();
        sd.metallic = m_metallic->eval(its).average();
        sd.specular = m_specular->eval(its).average();
        sd.specularTint = m_specularTint->eval(its).average();
        sd.roughness = m_roughness->eval(its).average();
        sd.anisotropic = m_anisotropic->eval(its).average();
        sd.sheen = m_sheen->eval(its).average();
        sd.sheenTint = m_sheenTint->eval(its).average();
        sd.clearcoat = m_clearcoat->eval(its).average();
        sd.clearcoatGloss = m_clearcoatGloss->eval(its).average();

        const auto lum = sd.baseColor.getLuminance();
        const Spectrum tintColor = lum > 0 ? (sd.baseColor / lum) : Spectrum(1);
        sd.specularColor = math::lerp(Spectrum(sd.metallic),
                                      Spectrum(sd.specular * 0.08f * math::lerp(Spectrum(sd.specularTint), Spectrum(1), tintColor)),
                                      sd.baseColor);

        return sd;
    }

    float diffuseProbability(const SurfaceData& sd, const BSDFSamplingRecord& bRec) const {
        const Vector wh = normalize((bRec.wi + bRec.wo));
        const auto cosThetaH = dot(bRec.wi, wh);
        const float fresnel = FrSchlick(sd.specularColor, cosThetaH).average();
        const float diffuse = (sd.baseColor * (1 - sd.metallic)).average() / M_PI;
        return diffuse / (diffuse + fresnel);
    }

    Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const override {
        bool hasSpecular = (bRec.typeMask & EGlossyReflection) &&
            (bRec.component == -1 || bRec.component == 0);
        bool hasDiffuse = (bRec.typeMask & EDiffuseReflection) &&
            (bRec.component == -1 || bRec.component == 1);

        if (measure != ESolidAngle ||
            Frame::cosTheta(bRec.wi) <= 0 ||
            Frame::cosTheta(bRec.wo) <= 0 ||
            (!hasSpecular && !hasDiffuse))
            return Spectrum(0.0f);

        const SurfaceData sd = prepareSurfaceData(bRec.its);

        Spectrum result(0.0f);

        // normalize lum. to isolate hue+sat
        const auto lum = sd.baseColor.getLuminance();
        const Spectrum tintColor = lum > 0 ? (sd.baseColor / lum) : Spectrum(1);
        const Spectrum specularColor = sd.specularColor;
        const Spectrum sheenColor = math::lerp(Spectrum(sd.sheenTint), Spectrum(1.0f), tintColor);

        if (hasDiffuse)
        {
            Spectrum diffuse = diffuseLobe(bRec, sd.roughness);
            Spectrum ss = sssLobe(bRec, sd.roughness);
            result += math::lerp(Spectrum(sd.subsurface), diffuse, ss) * sd.baseColor * (1.0f - sd.metallic);
        }

        if (hasSpecular)
        {
            result += sheenLobe(bRec, sd.sheen, sheenColor) * (1.0f - sd.metallic);
            result += specularLobe(bRec, sd.roughness, sd.anisotropic, specularColor);
            result += clearcoatLobe(bRec, sd.clearcoat, sd.clearcoatGloss);
        }

        return result * Frame::cosTheta(bRec.wo);
    }

    Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const override {
        if (!(bRec.typeMask & EDiffuseReflection)
            || !(bRec.typeMask & EGlossyReflection)
            || measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return 0.0f;

        const SurfaceData sd = prepareSurfaceData(bRec.its);

        const Vector3f wh = normalize(bRec.wi + bRec.wo);

        const auto diffuseProb = diffuseProbability(sd, bRec);
        const auto specularProb = 1 / (1 + sd.clearcoat);

        return diffuseProb * warp::squareToCosineHemispherePdf(bRec.wo) +
            (1 - diffuseProb) * specularProb * pdfGTR2Aniso(wh, sd.anisotropic, sd.roughness) / (4.0f * dot(wh, bRec.wo)) +
            (1 - diffuseProb) * (1 - specularProb) * pdfGTR1(wh, math::lerp(sd.clearcoatGloss, .1f, .001f)) / (4.0f * dot(wh, bRec.wo));
    }

    Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample_) const override {
        Float pdf;
        return sample(bRec, pdf, sample_);
    }

    Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf_, const Point2 &sample_) const override {
        if (!(bRec.typeMask & EDiffuseReflection) ||
            !(bRec.typeMask & EGlossyReflection) || 
            Frame::cosTheta(bRec.wi) <= 0)
            return Spectrum(0.0f);


        bRec.eta = 1.0f;
        bRec.sampledComponent = -1;
        bRec.sampledType = EDiffuseReflection | EGlossyReflection;

        const SurfaceData sd = prepareSurfaceData(bRec.its);

        bool sampleDiffuse = false;
        bool sampleSpecular = false;
        bool sampleClearcoat = false;

        Point2f sampleCopy = sample_;
        const auto diffuseProb = diffuseProbability(sd, bRec);
        const auto specularProb = 1 / (1 + sd.clearcoat);
        if (sampleCopy.y < diffuseProb)
        {
            sampleDiffuse = true;
            sampleCopy.y = math::clamp(sampleCopy.y / diffuseProb, 0.0f, 1.0f);
        }
        else
        {
            sampleCopy.y = math::clamp((sampleCopy.y - diffuseProb) / (1.0f - diffuseProb), 0.0f, 1.0f);

            if (sampleCopy.y < specularProb)
            {
                sampleSpecular = true;
                sampleCopy.y = math::clamp(sampleCopy.y / specularProb, 0.0f, 1.0f);
            }
            else
            {
                sampleClearcoat = true;
                sampleCopy.y = math::clamp((sampleCopy.y - specularProb) / (1.0f - specularProb), 0.0f, 1.0f);
            }
        }

        if (sampleDiffuse)
        {
            bRec.wo = warp::squareToCosineHemisphere(sampleCopy);
        }
        else
        {
            Vector3f wh;
            if (sampleSpecular)
                wh = sampleGTR2Aniso(sampleCopy, sd.anisotropic, sd.roughness);
            else if (sampleClearcoat)
                wh = sampleGTR1(sampleCopy, math::lerp(sd.clearcoatGloss, .1f, .001f));

            bRec.wo = normalize(reflect(bRec.wi, wh));
            if (Frame::cosTheta(bRec.wo) <= 0.0f)
            {
                return Spectrum{0.0f};
            }
        }

        pdf_ = pdf(bRec, EMeasure::ESolidAngle);
        return eval(bRec, EMeasure::ESolidAngle) / pdf_;
    }

    void addChild(const std::string &name, ConfigurableObject *child) override {
        if (child->getClass()->derivesFrom(MTS_CLASS(Texture))) {
            if (name == "baseColor")
                m_baseColor = static_cast<Texture *>(child);
            if (name == "subsurface")
                m_subsurface = static_cast<Texture *>(child);
            if (name == "metallic")
                m_metallic = static_cast<Texture *>(child);
            if (name == "specular")
                m_specular = static_cast<Texture *>(child);
            if (name == "specularTint")
                m_specularTint = static_cast<Texture *>(child);
            if (name == "roughness")
                m_roughness = static_cast<Texture *>(child);
            if (name == "anisotropic")
                m_anisotropic = static_cast<Texture *>(child);
            if (name == "sheen")
                m_sheen = static_cast<Texture *>(child);
            if (name == "sheenTint")
                m_sheenTint = static_cast<Texture *>(child);
            if (name == "clearcoat")
                m_clearcoat = static_cast<Texture *>(child);
            if (name == "clearcoatGloss")
                m_clearcoatGloss = static_cast<Texture *>(child);
        } else {
            BSDF::addChild(name, child);
        }
    }

    Float getRoughness(const Intersection &its, int component) const override {
        const Float roughness = m_roughness->eval(its).average();
        return roughnessToAlpha(roughness);
    }

    std::string toString() const override {
        std::ostringstream oss;
        oss << "Disney[" << endl
            << "  id = \"" << getID() << "\"," << endl
            << "  baseColor = " << indent(m_baseColor.toString()) << endl
            << "  subsurface = " << indent(m_subsurface.toString()) << endl
            << "  metallic = " << indent(m_metallic.toString()) << endl
            << "  specular = " << indent(m_specular.toString()) << endl
            << "  specularTint = " << indent(m_specularTint.toString()) << endl
            << "  roughness = " << indent(m_roughness.toString()) << endl
            << "  anisotropic = " << indent(m_anisotropic.toString()) << endl
            << "  sheen = " << indent(m_sheen.toString()) << endl
            << "  sheenTint = " << indent(m_sheenTint.toString()) << endl
            << "  clearcoat = " << indent(m_clearcoat.toString()) << endl
            << "  clearcoatGloss = " << indent(m_clearcoatGloss.toString()) << endl
            << "]";
        return oss.str();
    }

    Shader *createShader(Renderer *renderer) const override;

    MTS_DECLARE_CLASS()
private:
    ref<Texture> m_baseColor;
    ref<Texture> m_subsurface;
    ref<Texture> m_metallic;
    ref<Texture> m_specular;
    ref<Texture> m_specularTint;
    ref<Texture> m_roughness;           // perceptual linearly roughness
    ref<Texture> m_anisotropic;
    ref<Texture> m_sheen;
    ref<Texture> m_sheenTint;
    ref<Texture> m_clearcoat;
    ref<Texture> m_clearcoatGloss;
};

// ================ Hardware shader implementation ================

class DisneyShader : public Shader {
public:
    DisneyShader(Renderer *renderer, const Texture *reflectance)
        : Shader(renderer, EBSDFShader), m_reflectance(reflectance) {
        m_reflectanceShader = renderer->registerShaderForResource(m_reflectance.get());
    }

    bool isComplete() const {
        return m_reflectanceShader.get() != NULL;
    }

    void cleanup(Renderer *renderer) {
        renderer->unregisterShaderForResource(m_reflectance.get());
    }

    void putDependencies(std::vector<Shader *> &deps) {
        deps.push_back(m_reflectanceShader.get());
    }

    void generateCode(std::ostringstream &oss,
            const std::string &evalName,
            const std::vector<std::string> &depNames) const {
        oss << "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
            << "    if (cosTheta(wi) < 0.0 || cosTheta(wo) < 0.0)" << endl
            << "        return vec3(0.0);" << endl
            << "    return " << depNames[0] << "(uv) * inv_pi * cosTheta(wo);" << endl
            << "}" << endl
            << endl
            << "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
            << "    return " << evalName << "(uv, wi, wo);" << endl
            << "}" << endl;
    }

    MTS_DECLARE_CLASS()
private:
    ref<const Texture> m_reflectance;
    ref<Shader> m_reflectanceShader;
};

Shader *Disney::createShader(Renderer *renderer) const {
    return new DisneyShader(renderer, m_baseColor.get());
}

MTS_IMPLEMENT_CLASS(DisneyShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(Disney, false, BSDF)
MTS_EXPORT_PLUGIN(Disney, "Disney BRDF")
MTS_NAMESPACE_END
